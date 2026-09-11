# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:15:00 2026

@author: Benny


This module is used to generate synchronized multi-camera MP4 videos 
and GIF animations that follow individual particles or fibers along 
their 3D trajectories.

Can be run either through the workflow script (actions: 'trajectory_video'
or 'make_trajectory_video'), as a standalone script, or imported directly.
"""

from numpy import (array, zeros, arange, amax, amin, percentile,
                   degrees, arccos, arctan2, clip, linalg, hypot,
                   ceil, round as np_round, isnan)
from pandas import read_csv, DataFrame, merge
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
import os
import sys
import argparse
import imageio
from yaml import safe_load

# Try importing smart bbox decomposition if size_measure is installed
try:
    from geometry import decompose_blob_bbox, compute_local_scale
except ImportError:
    decompose_blob_bbox = None
    compute_local_scale = None




def resolve_filepath(path, base_dirs=None):
    '''
    Helper to resolve relative paths against candidate base directories.
    '''
    if path is None:
        return None
    if os.path.isabs(path) and os.path.exists(path):
        return path
        
    candidates = []
    if os.path.isabs(path):
        candidates.append(path)
    else:
        candidates.append(os.path.abspath(path))
        if base_dirs:
            for b in base_dirs:
                if b and os.path.exists(b):
                    candidates.append(os.path.abspath(os.path.join(b, path)))
                    
        for root in ['/Users/user/Desktop/Research',
                     '/Users/user/Desktop/Research/Data_Analysis/MyPTV_analysis']:
            candidates.append(os.path.abspath(os.path.join(root, path)))

    for cand in candidates:
        if os.path.exists(cand):
            return cand
    return candidates[0] if candidates else path




def load_frame_image(img_path):
    '''
    Loads an image from disk. Reads raw format (.dng) using rawpy
    or standard formats using imageio.
    '''
    ext = os.path.splitext(img_path)[1].lower()
    if ext == '.dng':
        import rawpy
        with rawpy.imread(img_path) as raw:
            return raw.raw_image.copy().astype(float)
    else:
        import imageio.v2 as iio
        im = iio.imread(img_path)
        if im.ndim == 3:
            im = 0.2989 * im[:,:,0] + 0.5870 * im[:,:,1] + 0.1140 * im[:,:,2]
        return im.astype(float)




class trajectory_video(object):
    '''
    A class used to generate a synchronized multi-camera video (MP4) 
    and animated GIF following a particle or fiber along its 3D trajectory.
    
    Can handle both spherical particles ('particles') and elongated fibers
    ('fibers'), with either the old bounding box method or the smart
    calibrated method.
    
    inputs -
    
    traj_file - string; path to the trajectory file (trajectories, 
                smoothed_trajectories, or .npz)
    camera_names - list or string; path to camera calibration files
    images_folder - string; directory containing camera image subfolders
    blob_files - list or string; paths to 2D blob coordinate files
    orientations_file - string (optional); path to fiber orientations file
    traj_id - int (optional); specific trajectory ID to follow
    shape - string; 'particles' or 'fibers'
    bbox_style - string; 'old' (direct dimensions) or 'smart' (decomposed)
    pad - int; crop half-width in pixels around particle
    fps - int; video framerate in frames per second
    fps_gif - int; GIF preview framerate
    save_mp4 - bool; if True, exports MP4 video
    save_gif - bool; if True, exports preview GIF
    out_mp4 - string (optional); output MP4 file path
    out_gif - string (optional); output GIF file path
    f_start, f_end - int (optional); frame range limits
    rec_name - string (optional); label for recording in figure title
    image_ext - string; extension of image files (default: '.dng')
    '''
    
    def __init__(self, traj_file=None, camera_names=None, images_folder=None,
                 blob_files=None, orientations_file=None, traj_id=None,
                 shape='particles', bbox_style='old', pad=40, fps=250,
                 fps_gif=10, save_mp4=True, save_gif=True, out_mp4=None,
                 out_gif=None, f_start=None, f_end=None, rec_name=None,
                 image_ext='.dng', base_dirs=None, trajectory_file=None,
                 orientation_file=None):
        
        self.traj_file = traj_file or trajectory_file
        self.camera_names = camera_names
        self.images_folder = images_folder
        self.blob_files = blob_files
        self.orientations_file = orientations_file or orientation_file
        self.traj_id = traj_id
        self.shape = str(shape).lower().strip()
        self.bbox_style = str(bbox_style).lower().strip()
        self.pad = int(pad)
        self.fps = int(fps)
        self.fps_gif = int(fps_gif)
        self.save_mp4 = bool(save_mp4)
        self.save_gif = bool(save_gif)
        self.out_mp4 = out_mp4
        self.out_gif = out_gif
        self.f_start = f_start
        self.f_end = f_end
        self.rec_name = rec_name
        self.image_ext = image_ext
        self.base_dirs = base_dirs or [os.getcwd()]

        # Validate trajectory file presence
        if not self.traj_file:
            raise ValueError('A trajectory_file must be specified.')

        # Validate dual file requirement for fibers
        if self.shape == 'fibers':
            is_npz = isinstance(self.traj_file, str) and self.traj_file.endswith('.npz')
            if not self.orientations_file and not is_npz:
                raise ValueError('For shape: fibers, both trajectory_file and orientation_file are required.')

        # 1) load cameras
        self.cams = self.load_cameras()
        self.num_cams = len(self.cams)

        # 2) locate image files
        self.cam_files, self.cam_dirs = self.load_images()

        # 3) load 2D blobs
        self.blob_dfs = self.load_blobs()

        # 4) load trajectory and orientations
        self.traj_data, self.traj_id_used, self.ori_map = self.load_trajectory()



    def load_cameras(self):
        '''
        Loads camera objects from calibration files.
        '''
        from myptv.imaging_mod import camera_wrapper
        
        if isinstance(self.camera_names, str):
            c_list = [c.strip() for c in self.camera_names.split(',') if c.strip()]
        else:
            c_list = list(self.camera_names)

        cams = []
        for cn in c_list:
            resolved = resolve_filepath(cn, self.base_dirs)
            if os.path.isdir(resolved):
                dir_path = resolved
                file_name = os.path.basename(cn) if os.path.basename(cn) != os.path.basename(resolved) else 'Cam1'
            else:
                dir_path, file_name = os.path.split(resolved)
                
            cam = camera_wrapper(file_name, dir_path)
            cam.load()
            cams.append(cam)
        return cams



    def load_images(self):
        '''
        Finds sorted image files for each camera folder.
        '''
        resolved = resolve_filepath(self.images_folder, self.base_dirs)
        if not os.path.exists(resolved):
            raise FileNotFoundError('Images folder not found: %s' % self.images_folder)

        base_folder = resolved
        last_elem = os.path.basename(os.path.normpath(resolved))
        if last_elem.lower().startswith('cam'):
            parent = os.path.dirname(os.path.normpath(resolved))
            if os.path.exists(parent):
                base_folder = parent

        cam_files = {}
        cam_dirs = {}
        for i in range(self.num_cams):
            candidate_names = ['Cam%d'%(i+1), 'cam%d'%(i+1), 'Camera_%d'%(i+1), 'camera_%d'%(i+1), 'Images_cam%d'%(i+1)]
            found_dir = None
            for cname in candidate_names:
                p = os.path.join(base_folder, cname)
                if os.path.isdir(p):
                    found_dir = p
                    break

            if found_dir is None:
                if self.num_cams == 1:
                    found_dir = base_folder
                else:
                    raise FileNotFoundError(f'Could not find folder for Camera {i+1} in {base_folder}. '
                                            f'Searching for folder names {candidate_names}')

            cam_dirs[i] = found_dir
            flist = sorted([
                f for f in os.listdir(found_dir)
                if f.lower().endswith(self.image_ext.lower()) and not f.startswith('._')
            ])
            if len(flist) == 0:
                raise FileNotFoundError('No %s images found in %s' % (self.image_ext, found_dir))
            cam_files[i] = flist

        return cam_files, cam_dirs



    def load_blobs(self):
        '''
        Loads 2D blob tables for all cameras into DataFrames.
        '''
        if self.blob_files is None:
            b_list = []
        elif isinstance(self.blob_files, str):
            b_list = [b.strip() for b in self.blob_files.split(',') if b.strip()]
        else:
            b_list = list(self.blob_files)

        blob_dfs = {}
        for i in range(self.num_cams):
            if i < len(b_list):
                bf = resolve_filepath(b_list[i], self.base_dirs)
                if os.path.exists(bf):
                    blob_dfs[i] = read_csv(bf, sep=r'\s+', header=None)
                    continue

            # Fallback search in working directories
            found = False
            for sdir in self.base_dirs:
                for pat in ['blobs_Cam%d_directions'%(i+1), 'blobs_cam%d_directions'%(i+1),
                            'blobs_Cam%d'%(i+1), 'blobs_cam%d'%(i+1)]:
                    cand = os.path.join(sdir, pat)
                    if os.path.exists(cand):
                        blob_dfs[i] = read_csv(cand, sep=r'\s+', header=None)
                        found = True
                        break
                if found:
                    break

            if not found:
                blob_dfs[i] = DataFrame()

        return blob_dfs



    def load_trajectory(self):
        '''
        Loads trajectory records and fiber orientation vectors.
        '''
        resolved_traj = resolve_filepath(self.traj_file, self.base_dirs)
        if not os.path.exists(resolved_traj):
            raise FileNotFoundError('Trajectory file not found: %s' % self.traj_file)

        orientations_map = {}
        if self.orientations_file:
            res_ori = resolve_filepath(self.orientations_file, self.base_dirs)
            if not os.path.exists(res_ori):
                raise FileNotFoundError('Orientation file not found: %s' % self.orientations_file)
            odf = read_csv(res_ori, sep=r'\s+', header=None)
            f_col = odf.columns[-1]
            for _, r in odf.iterrows():
                tid = int(r[0])
                frm = int(np_round(r[f_col]))
                orientations_map[(tid, frm)] = array([float(r[1]), float(r[2]), float(r[3])])

        # 1) Handle .npz format
        if resolved_traj.endswith('.npz'):
            import numpy as np
            npz_data = np.load(resolved_traj, allow_pickle=True)['data']
            tid = self.traj_id
            if tid is None:
                lengths = [len(tr) if tr is not None else 0 for tr in npz_data]
                tid = int(np.argmax(lengths))
            tr = npz_data[tid]
            return tr, tid, orientations_map

        # 2) Handle text trajectories file
        tdf = read_csv(resolved_traj, sep=r'\s+', header=None)
        valid_ids = [k for k in tdf[0].unique() if k > 0]
        if not valid_ids:
            valid_ids = list(tdf[0].unique())
        if not valid_ids:
            raise ValueError('No trajectory records found in %s' % resolved_traj)

        if self.traj_id is None or self.traj_id == 'longest':
            id_counts = tdf[tdf[0].isin(valid_ids)][0].value_counts()
            tid = int(id_counts.index[0])
        else:
            tid = int(self.traj_id)

        frame_col = tdf.columns[-1]
        sub_tr = tdf[tdf[0] == tid].sort_values(by=frame_col).copy()
        if len(sub_tr) == 0:
            raise ValueError('Trajectory ID %d not found in %s' % (tid, resolved_traj))
        sub_tr['frame_number'] = np_round(sub_tr[frame_col].values).astype(int)

        # Merge raw blob indices if smoothed trajectories file was given
        if len(tdf.columns) == 11 and 'smoothed' in os.path.basename(resolved_traj):
            raw_name = os.path.basename(resolved_traj).replace('smoothed_', '')
            raw_path = os.path.join(os.path.dirname(resolved_traj), raw_name)
            if os.path.exists(raw_path):
                try:
                    raw_df = read_csv(raw_path, sep=r'\s+', header=None)
                    raw_sub = raw_df[raw_df[0] == tid]
                    raw_f_col = raw_df.columns[-1]
                    cam_cols = [c for c in range(4, min(8, len(raw_df.columns) - 2))]
                    merged = merge(
                        sub_tr,
                        raw_sub[[raw_f_col] + cam_cols],
                        left_on='frame_number',
                        right_on=raw_f_col,
                        how='left',
                        suffixes=('', '_raw')
                    )
                    sub_tr = merged
                except Exception:
                    pass

        return sub_tr, tid, orientations_map



    def render(self):
        '''
        Renders the video frames and saves MP4 and GIF preview.
        '''
        import numpy as np
        
        if matplotlib.rcParams['text.usetex']:
            matplotlib.rcParams['text.usetex'] = False
        
        # Determine frame indices
        if isinstance(self.traj_data, np.ndarray):
            is_npz = True
            n_rows = len(self.traj_data)
            raw_frames = np_round(self.traj_data[:, 24] * 250.0).astype(int) if self.traj_data.shape[1] > 24 else arange(n_rows)
            frames = raw_frames % 1000000
        else:
            is_npz = False
            n_rows = len(self.traj_data)
            f_col = 'frame_number' if 'frame_number' in self.traj_data.columns else self.traj_data.columns[-1]
            frames = np_round(self.traj_data[f_col].values).astype(int)

        indices = []
        for idx in range(n_rows):
            frm = frames[idx]
            if self.f_start is not None and frm < self.f_start:
                continue
            if self.f_end is not None and frm > self.f_end:
                continue
            indices.append(idx)

        if not indices:
            raise ValueError('No frames to render for Trajectory %d in range [%s, %s]' % 
                             (self.traj_id_used, str(self.f_start), str(self.f_end)))

        start_frame = int(frames[indices[0]])
        end_frame = int(frames[indices[-1]])
        display_rec = self.rec_name or ('Traj %d' % self.traj_id_used)

        use_smart = (self.bbox_style == 'smart' and self.shape == 'fibers')
        if use_smart and decompose_blob_bbox is None:
            print('Notice: size_measure geometry module not found; using old bbox style.')
            use_smart = False
            self.bbox_style = 'old'

        print('--> Rendering Trajectory %d (%s, %s BBox, Frames %d..%d, %d frames)...' % 
              (self.traj_id_used, self.shape.upper(), self.bbox_style.upper(),
               start_frame, end_frame, len(indices)))

        video_frames = []

        if self.num_cams == 4:
            nrows, ncols = 2, 2
            figsize = (8.5, 8.5)
        elif self.num_cams == 2:
            nrows, ncols = 1, 2
            figsize = (9.0, 4.5)
        elif self.num_cams == 3:
            nrows, ncols = 1, 3
            figsize = (12.0, 4.2)
        else:
            nrows = int(ceil(self.num_cams / 2))
            ncols = 2
            figsize = (8.5, 4.2 * nrows)

        for row_idx in indices:
            f_idx = int(frames[row_idx])

            # 1) 3D coordinate and orientation
            if is_npz:
                pos_3d = self.traj_data[row_idx, 1:4]
                if self.traj_data.shape[1] > 12:
                    px, py, pz = self.traj_data[row_idx, 10:13]
                else:
                    px, py, pz = None, None, None
                blob_indices = [int(self.traj_data[row_idx, 19 + c_i]) if self.traj_data.shape[1] > 19 + c_i else -1
                                for c_i in range(self.num_cams)]
            else:
                row_vals = self.traj_data.iloc[row_idx].values
                pos_3d = array([float(row_vals[1]), float(row_vals[2]), float(row_vals[3])])
                
                if self.ori_map and (self.traj_id_used, f_idx) in self.ori_map:
                    px, py, pz = self.ori_map[(self.traj_id_used, f_idx)]
                elif self.shape == 'fibers' and len(row_vals) > 12:
                    px, py, pz = float(row_vals[10]), float(row_vals[11]), float(row_vals[12])
                else:
                    px, py, pz = None, None, None

                blob_indices = []
                for c_i in range(self.num_cams):
                    b_val = -1
                    if '%d_raw'%(4 + c_i) in self.traj_data.columns:
                        b_val = self.traj_data.iloc[row_idx]['%d_raw'%(4 + c_i)]
                    elif len(row_vals) >= 4 + self.num_cams and len(row_vals) != 11:
                        b_val = row_vals[4 + c_i]
                    blob_indices.append(int(b_val) if not isnan(b_val) else -1)

            has_ori = (px is not None and not isnan(px) and self.shape == 'fibers')
            if has_ori:
                phi_deg = degrees(arccos(clip(pz, -1.0, 1.0)))
                safe_px = px if abs(px) > 1e-12 else 1e-12
                theta_deg = degrees(arctan2(py, safe_px))
            else:
                phi_deg, theta_deg = 0.0, 0.0

            fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=100)
            fig.patch.set_facecolor('#1a1a1a')
            axes_flat = axes.flat if hasattr(axes, 'flat') else [axes]

            for c_i in range(self.num_cams):
                ax = axes_flat[c_i]
                ax.set_facecolor('#0a0a0a')

                flist = self.cam_files[c_i]
                img_file = flist[f_idx] if 0 <= f_idx < len(flist) else flist[min(max(0, f_idx), len(flist)-1)]
                img_path = os.path.join(self.cam_dirs[c_i], img_file)
                gray = load_frame_image(img_path)

                b_idx = blob_indices[c_i] if c_i < len(blob_indices) else -1
                matched = False
                c_row, c_col = None, None
                bw, bh = 6.0, 6.0
                mass = 0.0
                cos_th, sin_th = 1.0, 0.0

                df_blobs = self.blob_dfs.get(c_i, DataFrame())
                if not df_blobs.empty:
                    frame_blobs = df_blobs[df_blobs[5] == f_idx]
                    if b_idx >= 0 and b_idx < len(frame_blobs):
                        b_row = frame_blobs.iloc[b_idx]
                        c_row, c_col = float(b_row[0]), float(b_row[1])
                        bw, bh = float(b_row[2]), float(b_row[3])
                        mass = float(b_row[4])
                        if len(b_row) > 6:
                            cos_th = float(b_row[6])
                            sin_th = float(b_row[7])
                        matched = True
                    elif b_idx < 0 and not frame_blobs.empty:
                        proj = self.cams[c_i].projection(pos_3d)
                        p_col, p_row = float(proj[0]), float(proj[1])
                        dists = hypot(frame_blobs[1].values - p_col, frame_blobs[0].values - p_row)
                        min_idx = np.argmin(dists)
                        if dists[min_idx] <= 6.0:
                            b_row = frame_blobs.iloc[min_idx]
                            b_idx = min_idx
                            c_row, c_col = float(b_row[0]), float(b_row[1])
                            bw, bh = float(b_row[2]), float(b_row[3])
                            mass = float(b_row[4])
                            if len(b_row) > 6:
                                cos_th = float(b_row[6])
                                sin_th = float(b_row[7])
                            matched = True

                if not matched:
                    proj = self.cams[c_i].projection(pos_3d)
                    c_col, c_row = float(proj[0]), float(proj[1])

                # Crop around particle
                y_min = max(0, int(round(c_row - self.pad)))
                y_max = min(gray.shape[0], int(round(c_row + self.pad)))
                x_min = max(0, int(round(c_col - self.pad)))
                x_max = min(gray.shape[1], int(round(c_col + self.pad)))

                crop = gray[y_min:y_max, x_min:x_max]

                if crop.size > 0:
                    vmin, vmax = percentile(crop, [2, 99.8])
                    crop_norm = clip((crop - vmin) / (vmax - vmin + 1e-5), 0, 1)
                    ax.imshow(crop_norm, cmap='inferno', origin='upper')

                cx_loc = c_col - x_min
                cy_loc = c_row - y_min

                # Draw bounding box and markers
                if matched:
                    if self.shape == 'fibers':
                        u = array([sin_th, cos_th])
                        norm_u = linalg.norm(u)
                        u = u / norm_u if norm_u > 1e-6 else array([1.0, 0.0])
                        v = array([-u[1], u[0]])

                        if use_smart and decompose_blob_bbox is not None:
                            l_px, d_px = decompose_blob_bbox(bw, bh, cos_th, sin_th)
                            _, mm_per_px = compute_local_scale(self.cams[c_i], pos_3d)
                            l_mm = l_px * mm_per_px * 0.90
                            hL = (l_px * 0.90) / 2.0
                            hD = d_px / 2.0
                            status_str = 'Cam %d | L=%.1fpx (%.2fmm)' % (c_i+1, l_px, l_mm)
                        else:
                            hL = max(bw, bh) / 2.0
                            hD = min(bw, bh) / 2.0
                            status_str = 'Cam %d | Blob #%d' % (c_i+1, b_idx)

                        c1 = array([cx_loc, cy_loc]) + hL * u + hD * v
                        c2 = array([cx_loc, cy_loc]) + hL * u - hD * v
                        c3 = array([cx_loc, cy_loc]) - hL * u - hD * v
                        c4 = array([cx_loc, cy_loc]) - hL * u + hD * v

                        poly = Polygon([c1, c2, c3, c4], closed=True, edgecolor='cyan', facecolor='none', linewidth=1.8)
                        ax.add_patch(poly)
                        ax.plot(cx_loc, cy_loc, 'c+', markersize=6)
                        status_sub = 'Mass: %s' % f"{int(mass):,}"
                        status_color = '#00ffcc'

                    else:
                        h_w = max(bw / 2.0, 2.0)
                        h_h = max(bh / 2.0, 2.0)
                        p_corners = [
                            [cx_loc - h_w, cy_loc - h_h],
                            [cx_loc + h_w, cy_loc - h_h],
                            [cx_loc + h_w, cy_loc + h_h],
                            [cx_loc - h_w, cy_loc + h_h]
                        ]
                        poly = Polygon(p_corners, closed=True, edgecolor='cyan', facecolor='none', linewidth=1.8)
                        ax.add_patch(poly)
                        ax.plot(cx_loc, cy_loc, 'c+', markersize=6)
                        status_str = 'Cam %d | Blob #%d' % (c_i+1, b_idx)
                        status_sub = 'Mass: %s | %.0fx%.0fpx' % (f"{int(mass):,}", bw, bh)
                        status_color = '#00ffcc'

                else:
                    circle = Circle((cx_loc, cy_loc), 5, edgecolor='#ff4444', facecolor='none', linestyle='--', linewidth=1.5)
                    ax.add_patch(circle)
                    status_str = 'Cam %d | No Match' % (c_i+1)
                    status_sub = '3D Reproject'
                    status_color = '#ff6666'

                ax.set_xlim(0, 2 * self.pad)
                ax.set_ylim(2 * self.pad, 0)
                ax.axis('off')
                ax.set_title('%s\n%s' % (status_str, status_sub), color=status_color, fontsize=9, pad=3)

            for extra_ax in axes_flat[self.num_cams:]:
                extra_ax.axis('off')

            bbox_label = '(%s BBox)' % self.bbox_style.upper() if self.shape == 'fibers' else '(Particles)'
            header1 = '%s | Trajectory %d %s | Frame %03d / %d' % (display_rec, self.traj_id_used, bbox_label, f_idx, end_frame)
            if has_ori:
                header2 = ('3D: [%.1f, %.1f, %.1f] mm  |  p: [%+.2f, %+.2f, %+.2f]  |  θ: %+.1f°, φ: %.1f°' %
                           (pos_3d[0], pos_3d[1], pos_3d[2], px, py, pz, theta_deg, phi_deg))
            else:
                header2 = '3D: [%.2f, %.2f, %.2f] mm' % (pos_3d[0], pos_3d[1], pos_3d[2])

            fig.suptitle('%s\n%s' % (header1, header2), color='white', fontsize=11, weight='bold', y=0.98)
            plt.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.03, wspace=0.06, hspace=0.15)

            fig.canvas.draw()
            rgba = array(fig.canvas.buffer_rgba())[:, :, :3]
            video_frames.append(rgba)
            plt.close(fig)

        # 2) Save video and preview GIF
        if self.save_mp4 and self.out_mp4:
            os.makedirs(os.path.dirname(os.path.abspath(self.out_mp4)), exist_ok=True)
            print('Saving MP4 video at %d fps to %s...' % (self.fps, self.out_mp4))
            imageio.mimwrite(self.out_mp4, video_frames, fps=self.fps, codec='libx264', macro_block_size=None)
            print('Video saved: %s' % self.out_mp4)

        if self.save_gif and self.out_gif:
            os.makedirs(os.path.dirname(os.path.abspath(self.out_gif)), exist_ok=True)
            print('Saving preview GIF at %d fps to %s...' % (self.fps_gif, self.out_gif))
            imageio.mimwrite(self.out_gif, video_frames, fps=self.fps_gif, loop=0)
            print('GIF saved:   %s' % self.out_gif)

        return self.out_mp4, self.out_gif




def render_trajectory_video(traj_data, traj_id, cams, cam_files, cam_dirs,
                            blob_dfs, shape='particles', bbox_style='old',
                            orientations_map=None, pad=40, fps_mp4=250,
                            fps_gif=10, save_mp4=True, save_gif=True,
                            out_mp4=None, out_gif=None, frame_start=None,
                            frame_end=None, rec_name=None):
    '''
    Convenience function to render a trajectory video directly from arrays/objects.
    '''
    vid = trajectory_video.__new__(trajectory_video)
    vid.cams = cams
    vid.num_cams = len(cams)
    vid.cam_files = cam_files
    vid.cam_dirs = cam_dirs
    vid.blob_dfs = blob_dfs
    vid.shape = str(shape).lower().strip()
    vid.bbox_style = str(bbox_style).lower().strip()
    vid.pad = int(pad)
    vid.fps = int(fps_mp4)
    vid.fps_gif = int(fps_gif)
    vid.save_mp4 = bool(save_mp4)
    vid.save_gif = bool(save_gif)
    vid.out_mp4 = out_mp4
    vid.out_gif = out_gif
    vid.f_start = frame_start
    vid.f_end = frame_end
    vid.rec_name = rec_name
    vid.traj_data = traj_data
    vid.traj_id_used = traj_id
    vid.ori_map = orientations_map
    return vid.render()




def render_trajectory_video_from_params(params_file, **overrides):
    '''
    Parses a MyPTV YAML parameters file, extracts parameters with
    intelligent fallbacks across workflow blocks, and renders the video.
    '''
    params_path = os.path.abspath(params_file)
    if not os.path.exists(params_path):
        raise FileNotFoundError('Parameters file not found: %s' % params_file)

    with open(params_path, 'r') as f:
        yaml_data = safe_load(f)

    sections = {}
    if isinstance(yaml_data, list):
        for entry in yaml_data:
            if isinstance(entry, dict):
                for k, v in entry.items():
                    sections[k] = v if v is not None else {}
    elif isinstance(yaml_data, dict):
        sections = yaml_data

    vid_block = sections.get('trajectory_video') or sections.get('make_trajectory_video') or {}
    cfg = {**vid_block, **{k: v for k, v in overrides.items() if v is not None}}

    params_dir = os.path.dirname(params_path)
    base_dirs = [params_dir, os.getcwd(), os.path.dirname(params_dir)]

    # 1) Shape
    shape = cfg.get('shape')
    if not shape:
        seg_block = sections.get('segmentation', {})
        shape = seg_block.get('shape', 'particles')
    shape = str(shape).lower().strip()

    # 2) BBox style
    bbox_style = cfg.get('bbox_style', 'old')

    # 3) Camera calibration
    camera_names = cfg.get('camera_names')
    if not camera_names:
        camera_names = sections.get('matching', {}).get('camera_names')
    if not camera_names:
        camera_names = sections.get('analyze_calibration_error', {}).get('camera_names')
    if not camera_names:
        camera_names = sections.get('fiber_orientations', {}).get('camera_names')
    if not camera_names:
        camera_names = sections.get('calibration', {}).get('camera_name')
    if not camera_names:
        raise ValueError('No camera calibration (camera_names) found in parameters file.')

    # 4) Images folder and extension
    images_folder = cfg.get('images_folder')
    if not images_folder:
        seg_block = sections.get('segmentation', {})
        images_folder = seg_block.get('images_folder') or seg_block.get('images_folder1')
    if not images_folder:
        raise ValueError('No images_folder found in parameters file or arguments.')

    image_ext = cfg.get('image_extension') or sections.get('segmentation', {}).get('image_extension', '.dng')

    # 5) Trajectory file
    traj_file = cfg.get('trajectory_file') or cfg.get('traj_file')
    if not traj_file:
        traj_file = sections.get('smoothing', {}).get('save_name') or sections.get('smoothing', {}).get('trajectory_file')
    if not traj_file:
        traj_file = sections.get('tracking', {}).get('save_name')
    if not traj_file:
        traj_file = sections.get('stitching', {}).get('save_name')
    if not traj_file:
        raise ValueError('No trajectory file found in parameters file.')

    # 6) Blob files
    blob_files = cfg.get('blob_files')
    if not blob_files:
        if shape == 'fibers':
            blob_files = sections.get('fiber_orientations', {}).get('blob_files')
        if not blob_files:
            blob_files = sections.get('matching', {}).get('blob_files')

    # 7) Orientations file (fibers)
    orientations_file = cfg.get('orientation_file') or cfg.get('orientations_file')
    if not orientations_file and shape == 'fibers':
        orientations_file = sections.get('smoothed_orientations', {}).get('save_name') or sections.get('smoothed_orientations', {}).get('orientations_file')
        if not orientations_file:
            orientations_file = sections.get('fiber_orientations', {}).get('save_name')

    if shape == 'fibers':
        is_npz = isinstance(traj_file, str) and traj_file.endswith('.npz')
        if not orientations_file and not is_npz:
            raise ValueError('For shape: fibers, both trajectory_file and orientation_file are required.')

    # 8) Trajectory ID
    traj_id_req = cfg.get('traj_id') or cfg.get('traj_idx') or cfg.get('particle_id')

    # 9) Save paths
    save_dir = cfg.get('save_dir') or cfg.get('out_dir') or params_dir
    save_name = cfg.get('save_name')
    rec_name = cfg.get('rec_name')
    if not rec_name:
        for part in traj_file.replace('/', '_').split('_'):
            if part.lower().startswith('rec') and any(c.isdigit() for c in part):
                rec_name = part.capitalize()
                break
        if not rec_name:
            for part in images_folder.replace('/', '_').split('_'):
                if part.lower().startswith('rec') and any(c.isdigit() for c in part):
                    rec_name = part.capitalize()
                    break

    vid = trajectory_video(
        traj_file=traj_file,
        camera_names=camera_names,
        images_folder=images_folder,
        blob_files=blob_files,
        orientations_file=orientations_file,
        traj_id=traj_id_req,
        shape=shape,
        bbox_style=bbox_style,
        pad=int(cfg.get('pad', 40)),
        fps=int(cfg.get('fps', cfg.get('fps_mp4', 250))),
        fps_gif=int(cfg.get('fps_gif', 10)),
        save_mp4=bool(cfg.get('save_mp4', True)),
        save_gif=bool(cfg.get('save_gif', True)),
        f_start=cfg.get('frame_start') if cfg.get('frame_start') is None else int(cfg.get('frame_start')),
        f_end=cfg.get('frame_end') if cfg.get('frame_end') is None else int(cfg.get('frame_end')),
        rec_name=rec_name,
        image_ext=image_ext,
        base_dirs=base_dirs,
    )

    tag = ('_%s' % rec_name.lower()) if rec_name else ''
    default_prefix = 'traj%04d%s_%s_%s_bbox' % (vid.traj_id_used, tag, shape, bbox_style)
    if save_name:
        if save_name.endswith('.mp4') or save_name.endswith('.gif'):
            base_out = os.path.splitext(save_name)[0]
        else:
            base_out = os.path.join(save_dir, save_name)
    else:
        base_out = os.path.join(save_dir, default_prefix)

    vid.out_mp4 = '%s.mp4' % base_out
    vid.out_gif = '%s.gif' % base_out

    return vid.render()




def cli_main():
    '''
    Command-line interface for standalone execution.
    '''
    parser = argparse.ArgumentParser(
        description='Create synchronized multi-camera MP4 video & GIF for particle or fiber trajectory.'
    )
    parser.add_argument('params_file', help='Path to MyPTV parameters YAML file')
    parser.add_argument('--traj', '--traj-id', dest='traj_id', default=None, help="Trajectory ID to render (or 'longest')")
    parser.add_argument('--shape', choices=['particles', 'fibers'], default=None, help="Target shape ('particles' or 'fibers')")
    parser.add_argument('--bbox', '--bbox-style', dest='bbox_style', choices=['old', 'smart'], default=None, help='Bounding box style')
    parser.add_argument('--pad', type=int, default=None, help='Crop region half-size in pixels (default: 40)')
    parser.add_argument('--fps', type=int, default=None, help='MP4 framerate (default: 250)')
    parser.add_argument('--fps-gif', type=int, default=None, help='GIF framerate (default: 10)')
    parser.add_argument('--no-gif', action='store_true', help='Disable GIF preview export')
    parser.add_argument('--images-folder', default=None, help='Path to images directory (overrides params file)')
    parser.add_argument('--out-dir', default=None, help='Output directory for MP4 and GIF')
    parser.add_argument('--save-name', default=None, help='Custom output filename prefix')
    parser.add_argument('--traj-file', '--trajectory-file', dest='traj_file', default=None, help='Path to trajectory file')
    parser.add_argument('--orientation-file', '--orientations-file', dest='orientation_file', default=None, help='Path to orientation file (fibers)')
    parser.add_argument('--blob-files', default=None, help='Paths to 2D blob files')
    parser.add_argument('--f-start', type=int, default=None, help='Starting frame number')
    parser.add_argument('--f-end', type=int, default=None, help='Ending frame number')

    args = parser.parse_args()

    render_trajectory_video_from_params(
        args.params_file,
        trajectory_file=args.traj_file,
        orientation_file=args.orientation_file,
        blob_files=args.blob_files,
        traj_id=args.traj_id,
        shape=args.shape,
        bbox_style=args.bbox_style,
        pad=args.pad,
        images_folder=args.images_folder,
        fps=args.fps,
        fps_gif=args.fps_gif,
        save_gif=(not args.no_gif),
        save_dir=args.out_dir,
        save_name=args.save_name,
        frame_start=args.f_start,
        frame_end=args.f_end,
    )


if __name__ == '__main__':
    cli_main()
