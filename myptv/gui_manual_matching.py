# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A Gui that is intended to help users perform manual
matching by clicking points over images. This can be
used to analyze the quality of the calibration or to
easily extract the 3D location of certain features 
in images.  
"""

from myptv.imaging_mod import camera
from myptv.segmentation_mod import particle_segmentation
from myptv.calibrate_mod import calibrate
from myptv.utils import match_calibration_blobs_and_points
from PIL import Image, ImageTk
from tkinter import Label,Canvas,LabelFrame,Entry,Tk,Scrollbar,Button,Listbox,END
from numpy import array
import os






class initial_cal_gui(object):
    '''
    This is a Tkinter based graphical user interface that can be used to mark
    points on a static image, give their lab space coordinates, and then save
    the data as a text file in the format used by the calibration processes.
    '''
    
    
    def __init__(self, cam_list, image_list, cameras_folder = '.'):
        '''
        input: 
        
        cam_list - list of string, camera names
        
        image_list - list of strings, path to the images to be used. The
                     images should be ordered the same as the cameras in 
                     cam_list
                     
        cameras_folder - string, path of the folder that holds the camera 
                         files. default is in the current working directory.
        
        '''
        self.image_names = image_list
        self.images = [Image.open(img) for img in self.image_names]
        self.images_size = [img.size for img in self.images]
        
        
        self.cam_names = cam_list
        self.cam_res = [image_size for image_size in self.images_size]
        self.cam_list = []
        
        for e,cn in enumerate(self.cam_names):
            cam = camera(cn, self.cam_res[e])
            cam.load(cameras_folder)
            self.cam_list.append(cam)
        
        
        # try to find a calibration folder
        ls = os.listdir('.')
        self.folder = '.'
        for fname in ls:
            if fname in ['calibration', 'Calibration', 'cal', 'Cal']:
                if os.path.isdir(os.path.join('.', fname)):
                    self.folder = os.path.join('.', fname)
        
        
        self.xy_marked = (-1, -1)
        self.point_list = []      # <-- list of points to save
        self.point_markers = []   # <-- position of crosses
        self.z = 1.0              # <-- Zoom level
        
        # set the window
        self.root = Tk()
        self.root.geometry('1100x700+320+70')
        self.root.title('MyPTV: Initial Calibration GUI')
        #self.root.configure(background='#709eba')

        
        # =====================================================================
        # (1) Image frame, first column
        
        # load the image
        photo = ImageTk.PhotoImage(self.images[0])
        
        
        
        # place the image inside a canvas in a frame
        self.board = Label(image=photo)
        frame = LabelFrame(self.root, bg='#c2c2c2', width=700, height=450,
                           padx=(5), pady=5)
        frame.grid(row=0, column=0, padx=(5), pady=5, sticky='nsew')
        #frame.place(x=120, y=30)
        
        sr = (0,0,self.images_size[0][0],self.images_size[0][1])
        self.canvas = Canvas(frame, scrollregion=sr, height=450, width=700)       
        self.canvas.grid(row=0,column=1, sticky='ewns', padx=(5), pady=5)
        #self.canvas.focus()
        
        # add scrollbars to the image
        self.hbar=Scrollbar(frame,orient='horizontal')
        self.hbar.grid(row=1,column=1, sticky='ew')
        self.hbar.config(command=self.canvas.xview)
        self.vbar=Scrollbar(frame,orient='vertical')
        self.vbar.grid(row=0,column=0, sticky='nse')
        self.vbar.config(command=self.canvas.yview)
        
        self.canvas.config(xscrollcommand=self.hbar.set,
                            yscrollcommand=self.vbar.set)
        self.canvas.create_image(0, 0, image=photo, anchor='nw')
        
        
        
        # setup click bottun - when we click on the image it take the position
        self.canvas.bind("<Button-1>", self.location_handler)
        
        # zoom in and zoom out by pressing + and -
        self.root.bind('+', self.zoomIn)
        self.root.bind('-', self.zoomOut)
        
        # move cross with the arrow keys
        self.root.bind('<Shift-Left>', self.leftKey)
        self.root.bind('<Shift-Right>', self.rightKey)
        self.root.bind('<Shift-Up>', self.upKey)
        self.root.bind('<Shift-Down>', self.downKey)
        
        
        
        
        
        
        
        
        
        # ====================================================================
        # (3) second column
        
        
        Column3 = LabelFrame(self.root, padx=2, pady=10, width=100, 
                                  bg='#c2c2c2')
        Column3.grid(row=0, column=3, padx=(2), pady=10, sticky='nsew')
        
        
        #==================================
        # Marking points for initial calibration
        
        # Buttons frame
        button_frame = LabelFrame(Column3, text='3) mark image points', 
                                  padx=2, pady=8, width=100)
        button_frame.grid(row=0, column=0, columnspan=2, sticky='nwe', padx=2, 
                          pady=8)
        
        # add point button
        add_button = Button(button_frame, text='Mark point', 
                                command = self.addPoint, padx=2, pady=4)
        add_button.grid(row=4, column=0, padx=2, pady=2, sticky='ew')
        
        # forget last point button
        forget_last_button = Button(button_frame, text='Forget last point', 
                                command = self.forgetLast, padx=2, pady=4)
        forget_last_button.grid(row=5, column=0, padx=2, pady=2, sticky='ew')
        
        
        
        
        # camera selector
        camera_select_frame = LabelFrame(button_frame, padx=2, pady=2)
        camera_select_frame.grid(row=0, column=0, columnspan=2, sticky='ew',
                          padx=2, pady=2)
        self.listbox = Listbox(camera_select_frame)
        self.listbox.grid(row=0, column=0, rowspan=1, sticky='w', padx=2, pady=2)
        for cn in self.cam_names:
            self.listbox.insert(END, cn)
        
        
        
        # selection indicators
        select_frame = LabelFrame(button_frame, padx=2, pady=2)
        select_frame.grid(row=1, column=0, columnspan=2, sticky='ew',
                          padx=2, pady=2)
        self.xloc = Label(select_frame, text='x image:', padx=2, pady=2)
        self.yloc = Label(select_frame, text='y image:', padx=2, pady=2)
        self.Xloc = Label(select_frame, text='-', padx=2, pady=2)
        self.Yloc = Label(select_frame, text='-', padx=2, pady=2)
        self.xloc.grid(row=2, column=0, rowspan=1, sticky='w', padx=2, pady=2)
        self.yloc.grid(row=3, column=0, rowspan=1, sticky='w', padx=2, pady=2)
        self.Xloc.grid(row=2, column=1, sticky='w', padx=2, pady=2)
        self.Yloc.grid(row=3, column=1, sticky='w', padx=2, pady=2)
        
        
        # lab space coordinates textboxes
        lab_frame = LabelFrame(button_frame, padx=2, pady=2)
        lab_frame.grid(row=8, column=0, columnspan=2, sticky='ew',
                       padx=2, pady=2)
        self.x_input_label = Label(lab_frame, text='x lab:', padx=2, pady=2)
        self.y_input_label = Label(lab_frame, text='y lab:', padx=2, pady=2)
        self.z_input_label = Label(lab_frame, text='z lab:', padx=2, pady=2)
        self.x_input = Entry(lab_frame, width=9)
        self.x_input.insert(0,'0.0')
        self.y_input = Entry(lab_frame, width=9)
        self.y_input.insert(0,'0.0')
        self.z_input = Entry(lab_frame, width=9)
        self.z_input.insert(0,'0.0')
        
        self.x_input_label.grid(row=4, column=0, sticky='w', padx=2, pady=2)
        self.y_input_label.grid(row=5, column=0, sticky='w', padx=2, pady=2)
        self.z_input_label.grid(row=6, column=0, sticky='w', padx=2, pady=2)
        self.x_input.grid(row=4, column=1, sticky='w', padx=2, pady=2)
        self.y_input.grid(row=5, column=1, sticky='w', padx=2, pady=2)
        self.z_input.grid(row=6, column=1, sticky='w', padx=2, pady=2)
        
        
        # set the mouse motion printing
        self.canvas.bind("<Motion>", self.motion)
        Mousepos_frame = LabelFrame(button_frame, padx=2, pady=2)
        Mousepos_frame.grid(row=3, sticky='we', columnspan=2, padx=2, pady=2)
        self.Mousepos = Label(Mousepos_frame, text=' - , - ', padx=2, pady=2)
        self.Mousepos.grid(row=0, column=0, sticky='e', padx=2, pady=2)
        
        
        
        
        # quit button
        
        quit_frame = LabelFrame(Column3, text='', 
                                  padx=2, pady=8, width=100)
        quit_frame.grid(row=2, column=0, columnspan=2, sticky='sew', padx=2, 
                          pady=2)
        
        quit_button = Button(quit_frame, text='Quit', width=19,
                                command = self.Quit, padx=2, pady=2) 
        quit_button.grid(row=0, column=0, padx=2, pady=2, sticky='sew')
        
        
        
        #=====================================================================
        # configure hte frames and run main loop
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        Column3.rowconfigure(0, weight=1)
        Column3.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        self.root.mainloop()
        
        
        
        
    
    def calibrate(self):
        '''Calibrates the camera using marked points'''
        
        cpf = os.path.join(self.folder, self.cam_name+'_manualPoints')
        self.cam = camera(self.cam_name, self.cam_res, cal_points_fname = cpf)
        self.cam.load('.')
        print('camera data loaded successfully.')
        cal = calibrate(self.cam, self.cam.lab_points, self.cam.image_points)
        err = cal.mean_squared_err()
        print('initial error: %.3f pixels\n'%(err))
        
        print('\nAttempting to minimize calibration error...\n')
        err_iminus1 = err
        for i in range(10):
            cal.searchCalibration()
            err = cal.mean_squared_err()
            
            print(i, err)
            
            if err<0.5:
                break
            
            if abs(err - err_iminus1)/err_iminus1 < 0.001:
                break
        
            err_iminus1 = err
        
        self.error_input.config(text = '%0.3f'%err)
        self.cam.save('.')
        print('Finished with error: %.3f pixels\n'%(err))
            
            
        
        
    def show_calibration(self):
        '''Plots the calibration points images on the calibration image'''
        
        # refresh the image
        image = Image.open(self.image_name)
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.ANTIALIAS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        
        # plot point projections
        for p in self.cam.lab_points:
            proj = self.cam.projection(p)
            x0 = int(int(proj[0])*self.z) - int(self.hbar.get()[1])
            y0 = int(int(proj[1])*self.z) - int(self.vbar.get()[1])
            
            print(proj)
            self.canvas.create_oval(x0-4, y0-4, x0+4, y0+4, 
                                    fill='#cf9417')
        
    
    
    
    def saveTargetPoints(self):
        '''Will save the matched segmented blobs and targets, making a 
           calibration points file.'''
        saveName = os.path.join(self.folder, self.cam_name + '_cal_points')
        self.mtf.save_results(saveName)
        print('Calibration point file "%s" saved'%saveName)
        
        
        
        
    def matchTargetPoints(self):
        '''Will match the target points to the segmented blobs'''
        
        segmented_file = os.path.join(self.folder, self.cam_name + '_CalBlobs')
        self.mtf = match_calibration_blobs_and_points(self.cam, segmented_file,
                                                      self.target_fname)
        self.mtf.pair_points()
    
        # plot the pairs
        self.plot_matched_point_pair()
      
        
      
            
    def plot_matched_point_pair(self):
        '''
        Will plot the pairs of calibration points (segmented and their 
        projections).
        '''
        self.segmented = []
        for pair in self.mtf.point_pairs:
            # the segmented image space calibration points
            xImg = int(pair[0]*self.z) - int(self.hbar.get()[1]) 
            yImg = int(pair[1]*self.z) - int(self.vbar.get()[1])
            
            wx = 4*self.z; wy = 4*self.z
            self.canvas.create_rectangle(xImg-wx, yImg-wy, xImg+wx, yImg+wy, 
                                         outline="#2bcf0e", width=1)
            
            # the projection of the lab space points 
            eta, zeta = self.cam.projection(pair[2:])
            xProj = int(eta*self.z) - int(self.hbar.get()[1])
            yProj = int(zeta*self.z) - int(self.vbar.get()[1])
            
            wx = 2*self.z; wy = 2*self.z
            self.canvas.create_oval(xProj-wx, yProj-wy, xProj+wx, yProj+wy, 
                                         fill="#d94559", width=0)
            
            # connect the pairs
            self.canvas.create_line(xImg, yImg, xProj, yProj, 
                                    fill="#991527", width=1)
            
            


    def genCamFile(self):
        '''
        Generates a camera file with the given initial guess
        '''
        ox = float(self.Ox_input.get()) 
        oy = float(self.Oy_input.get()) 
        oz = float(self.Oz_input.get())
        tx = float(self.xori_input.get()) 
        ty = float(self.yori_input.get()) 
        tz = float(self.zori_input.get()) 
        f = float(self.f_input.get()) 
        #xh = float(self.xh_input.get()) 
        #yh = float(self.yh_input.get()) 
        
        self.cam.O = [ox,oy,oz]
        self.cam.theta = [tx,ty,tz]
        self.cam.f = f
        #self.cam.xh = xh
        #self.cam.yh = yh
        self.cam.calc_R()
        self.cam.save('.')
        
        print(self.cam)
        print('\nCamera file generated. \n')
        return None




    def sementImage(self):
        '''Segments the image and save the blob file'''
        from numpy import zeros
        
        x0,y0 = int(self.ROIx0.get()), int(self.ROIy0.get())
        x1,y1 = int(self.ROIx1.get()), int(self.ROIy1.get())
        mask = zeros((self.cam_res[1], self.cam_res[0]))
        mask[y0:y1,x0:x1] = 1
        
        self.segmentationParams ={'image': array(self.image),
                                  'threshold': float(self.threshold_input.get()), 
                                  'max_xsize': float(self.xmax_input.get()), 
                                  'max_ysize': float(self.ymax_input.get()), 
                                  'min_xsize': float(self.xmin_input.get()), 
                                  'min_ysize': float(self.ymin_input.get()),
                                  'min_mass': float(self.minMass_input.get()),
                                  'max_mass': float(self.maxMass_input.get()),
                                  'sigma': float(self.sigma_input.get()), 
                                  'median': int(self.median_input.get()),
                                  'local_filter': int(self.local_input.get()), 
                                  'mask': mask,
                                  'method': 'labeling',
                                  'particle_size':8}
        
        for k in ['sigma', 'median', 'local_filter']:
            if self.segmentationParams[k] == 0.0:
                self.segmentationParams[k] = None
        
        self.particleSegment = particle_segmentation(**self.segmentationParams)
        self.particleSegment.get_blobs()
        self.particleSegment.apply_blobs_size_filter()
        
        print('\nSegmenting image...\n')
        print('blobs found:', len(self.particleSegment.blobs))
        
        self.segmented = []
        for b in self.particleSegment.blobs:
            self.segmented.append((b[0][1], b[0][0], b[1][1], b[1][0]))
            
        # plot the segmented particles over a refreshed image
        image = Image.open(self.image_name)
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.ANTIALIAS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))

        self.plotSegmented()
        
        return None
    
    
    def save_blobs(self):
        '''
        Saves the results of the segmentation
        '''
        saveName = os.path.join(self.folder, self.cam_name + '_CalBlobs')
        self.particleSegment.save_results(saveName)
        print('\nFile saved. Done.\n')
    
    
    def plotSegmented(self):
        for x_, y_, wx, wy in self.segmented:
            x_ = int(x_*self.z) - int(self.hbar.get()[1])
            y_ = int(y_*self.z) - int(self.vbar.get()[1])
            wx = wx/2*self.z; wy = wy/2*self.z
            self.canvas.create_rectangle(x_-wx, y_-wy, x_+wx, y_+wy, 
                                         outline="#2bcf0e", width=1)
    

    def addPoint(self):
        '''add marked point to list'''
        x_im, y_im = self.xy_marked
        x_lab = float(self.x_input.get())
        y_lab = float(self.y_input.get())
        z_lab = float(self.z_input.get())
        self.point_list.append([x_im, y_im, x_lab, y_lab, z_lab])
        
        print('Added to list: ', self.point_list[-1])
        self.xy_marked = (-1, -1)
        self.mark_points()
            
        
    def forgetLast(self):
        p = self.point_list.pop(-1)
        print('Deleted: ', p)
        self.mark_points()
        
        
    def rightKey(self, event):
        '''right key = move cross'''
        self.xy_marked = (self.xy_marked[0]+1, self.xy_marked[1])
        self.mark_points()
        self.Xloc.configure(text = self.xy_marked[0]) 
        
        
    def leftKey(self, event):
        '''left key = move cross'''
        self.xy_marked = (self.xy_marked[0]-1, self.xy_marked[1])
        self.mark_points()
        self.Xloc.configure(text = self.xy_marked[0])
        
    
    def downKey(self, event):
        '''left key = move cross'''
        self.xy_marked = (self.xy_marked[0], self.xy_marked[1]+1)
        self.mark_points()
        self.Yloc.configure(text = self.xy_marked[1])
        
        
    def upKey(self, event):
        '''left key = move cross'''
        self.xy_marked = (self.xy_marked[0], self.xy_marked[1]-1)
        self.mark_points()
        self.Yloc.configure(text = self.xy_marked[1])
        
    
    
    def zoomIn(self, event):
        '''zoom in the image with + key'''
        self.z = self.z*1.15
        image = Image.open(self.image_name)
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.ANTIALIAS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        self.mark_points()
        
        if hasattr(self, 'mtf'):
            self.plot_matched_point_pair()
        
        else:
            self.plotSegmented()
            
    
        
    def zoomOut(self, event):
        '''zoom out with - key'''
        self.z = self.z*(1/1.15)
        image = Image.open(self.image_name)
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.ANTIALIAS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        self.mark_points()
        
        if hasattr(self, 'mtf'):
            self.plot_matched_point_pair()
        
        else:
            self.plotSegmented()
        
        
    def location_handler(self, event):
        '''handling the location of the mouse pointer'''
        #x,y = int(event.x/z), int(event.y/z)
        
        x = int(self.canvas.canvasx(event.x)/self.z) + int(self.hbar.get()[1])
        y = int(self.canvas.canvasy(event.y)/self.z) + int(self.vbar.get()[1])
        
        self.Xloc.configure(text = x) 
        self.Yloc.configure(text = y)
        self.xy_marked = (x, y)
        #tracking_dic[i%len(imagelist)]=(x,y)
        self.mark_points()
        print( x,y )
        
    
    def motion(self, event):
        '''handling the location of the mouse pointer'''
        #x, y = int(event.x/z), int(event.y/z)
        x = int(self.canvas.canvasx(event.x)/self.z) + int(self.hbar.get()[1])
        y = int(self.canvas.canvasy(event.y)/self.z) + int(self.vbar.get()[1])
        self.Mousepos.configure(text = '(%d, %d)'%(x, y))   
        
    
    def mark_points(self):
        '''Puts crosses on the location of saves trajectory points'''
        
        # delete the old cross
        for c in self.point_markers:
            self.canvas.delete(c)
        
        x, y = self.xy_marked
        x_ = int(x*self.z) - int(self.hbar.get()[1])
        y_ = int(y*self.z) - int(self.vbar.get()[1])
        c1 = self.canvas.create_line(x_-4, y_, x_+4, y_, 
                                     fill="#e31010", width=1)
        c2 = self.canvas.create_line(x_, y_-4, x_, y_+4, 
                                     fill="#e31010", width=1)    
        self.point_markers.append(c1)
        self.point_markers.append(c2)
        
        for p in self.point_list:
            x, y = p[0], p[1]
            x_ = int(x*self.z) - int(self.hbar.get()[1])
            y_ = int(y*self.z) - int(self.vbar.get()[1])
            c1 = self.canvas.create_line(x_-4, y_, x_+4, y_, 
                                         fill="#001ced", width=1)
            c2 = self.canvas.create_line(x_, y_-4, x_, y_+4, 
                                         fill="#001ced", width=1)    
            self.point_markers.append(c1)
            self.point_markers.append(c2)
        
    
    def Save(self):
        '''save the manually selected calibration points'''
        from numpy import savetxt
        saveName = os.path.join(self.folder, self.cam_name+'_manualPoints')
        savetxt(saveName, self.point_list, 
                fmt='%.1f', delimiter='\t')
        print('Points saved at: %s'%saveName)
        
        
    def Quit(self):
        '''quit the app'''
        self.root.destroy()




if __name__ == '__main__':
    im_fname = ['../example/Calibration/cal1.tif',
                '../example/Calibration/cal2.tif',
                '../example/Calibration/cal3.tif']
    
    cameras_folder = '../example'
    
    camera_names = ['cam1', 'cam2', 'cam3']
    
    gui = initial_cal_gui(camera_names, im_fname, cameras_folder=cameras_folder)


