# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A gui for the initial calibration of a Tsai-model camera. 
"""

from myptv.TsaiModel.camera import camera_Tsai
from myptv.TsaiModel.calibrate import calibrate_Tsai
from myptv.segmentation_mod import particle_segmentation
from myptv.utils import match_calibration_blobs_and_points

from PIL import Image, ImageTk
from tkinter import Label, Canvas, LabelFrame, Entry, Tk, Scrollbar, Button
from numpy import array
import os






class initial_cal_gui(object):
    '''
    This is a Tkinter based graphical user interface that can be used to mark
    points on a static image, give their lab space coordinates, and then save
    the data as a text file in the format used by the calibration processes.
    '''
    
    
    def __init__(self, cam_name, image_name, target_fname):
        '''
        input: 
        
        cam_name - string, the name of the cameras that we wish to calibrate.
        
        image_name - the calibration image used choose the points.
        
        target_fname - the file name of the calibration target's target file.
        '''
        self.image_name = image_name
        self.image = Image.open(self.image_name)
        image_size = self.image.size
        
        self.target_fname = target_fname
        
        self.cam_name = cam_name
        self.cam_res = image_size
        try:
            self.cam = camera_Tsai(self.cam_name, 
                          cal_points_fname = self.cam_name+'_manualPoints')
            self.cam.load('.')
        except:
            self.cam = camera_Tsai(self.cam_name)

        self.segmented = []
        
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
        photo = ImageTk.PhotoImage(self.image)
        
        
        
        # place the image inside a canvas in a frame
        self.board = Label(image=photo)
        frame = LabelFrame(self.root, bg='#c2c2c2', width=700, height=450,
                           padx=(5), pady=5)
        frame.grid(row=0, column=0, padx=(5), pady=5, sticky='nsew')
        #frame.place(x=120, y=30)
        
        sr = (0,0,image_size[0],image_size[1])
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
        # (2) The second column
        second_column = LabelFrame(self.root, padx=2, pady=10, width=100, 
                                   bg='#c2c2c2')
        second_column.grid(row=0, column=1, padx=(2), pady=10, sticky='nsew')
        
        
        
        # =========================
        # cam file generation frame
        camFile_frame = LabelFrame(second_column, padx=2, pady=8, width=100,
                                   text='1) cam file generation')
        camFile_frame.grid(row=0, column=0, padx=5, pady=8, sticky='nwe')
        
        
        
        name_dashboard = LabelFrame(camFile_frame)
        name_dashboard.grid(row=0, column=0, columnspan=1, sticky='sw', padx=5, 
                          pady=5)
        
        self.camName = Label(name_dashboard, text='Camera:', padx=2, pady=2)
        self.camName.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.camName_input = Label(name_dashboard, text=self.cam_name, padx=2, pady=2,
                                 width=10, bg='white')
        self.camName_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        pos_dashboard = LabelFrame(camFile_frame)
        pos_dashboard.grid(row=1, column=0, columnspan=1, sticky='nsw', padx=5, 
                          pady=2)
        
        self.xloc = Label(pos_dashboard, text='Ox:', padx=2, pady=2)
        self.xloc.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.Ox_input = Entry(pos_dashboard, width=14)
        self.Ox_input.insert(0,'0.0')
        self.Ox_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.yloc = Label(pos_dashboard, text='Oy:', padx=2, pady=2)
        self.yloc.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.Oy_input = Entry(pos_dashboard, width=14)
        self.Oy_input.insert(0,'0.0')
        self.Oy_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.zloc = Label(pos_dashboard, text='Oz:', padx=2, pady=2)
        self.zloc.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.Oz_input = Entry(pos_dashboard, width=14)
        self.Oz_input.insert(0,'0.0')
        self.Oz_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        
        # ori_Dashboard frame - where the camera angles are shown
        ori_dashboard = LabelFrame(camFile_frame)
        ori_dashboard.grid(row=2, column=0, columnspan=1, sticky='nw', padx=5, 
                          pady=2)
        
        self.xori = Label(ori_dashboard, text=chr(1012)+'x:', padx=2, pady=2)
        self.xori.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.xori_input = Entry(ori_dashboard, width=14)
        self.xori_input.insert(0,'0.0')
        self.xori_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.yori = Label(ori_dashboard, text=chr(1012)+'y:', padx=2, pady=2)
        self.yori.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.yori_input = Entry(ori_dashboard, width=14)
        self.yori_input.insert(0,'0.0')
        self.yori_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.zori = Label(ori_dashboard, text=chr(1012)+'z:', padx=2, pady=2)
        self.zori.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.zori_input = Entry(ori_dashboard, width=14)
        self.zori_input.insert(0,'0.0')
        self.zori_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        # int_dashboard frame - where focus and xh, yh are shown
        int_dashboard = LabelFrame(camFile_frame)
        int_dashboard.grid(row=3, column=0, columnspan=1, sticky='nw', padx=5, 
                          pady=2)
        
        self.f = Label(int_dashboard, text='f:', padx=2, pady=2)
        self.f.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.f_input = Entry(int_dashboard, width=14)
        self.f_input.insert(0,'0.0')
        self.f_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        # self.xh = Label(int_dashboard, text='xh:', padx=2, pady=2)
        # self.xh.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        # self.xh_input = Entry(int_dashboard, width=14)
        # self.xh_input.insert(0,'0.0')
        # self.xh_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        # self.yh = Label(int_dashboard, text='yh:', padx=2, pady=2)
        # self.yh.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        # self.yh_input = Entry(int_dashboard, width=14)
        # self.yh_input.insert(0,'0.0')
        # self.yh_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        gen_cam_button = Button(camFile_frame, text='Generate cam file', 
                                command = self.genCamFile, padx=2, pady=7)
        gen_cam_button.grid(row=5, column=0, padx=2, pady=7, sticky='ew')
        
        
        
        
        
        
        
        
        # ==================
        # segmentation frame
        segmentation_frame = LabelFrame(second_column, padx=2, pady=8, 
                                        width=100, text='2) image segmentation')
        segmentation_frame.grid(row=1, column=0, padx=5, pady=8, sticky='s')
        
        
        
        self.threshold = Label(segmentation_frame, text='threshold:', padx=2, pady=2)
        self.threshold.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.threshold_input = Entry(segmentation_frame, width=7)
        self.threshold_input.insert(0,'20')
        self.threshold_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        self.xmin = Label(segmentation_frame, text='x min size:', padx=2, pady=2)
        self.xmin.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.xmin_input = Entry(segmentation_frame, width=7)
        self.xmin_input.insert(0,'2')
        self.xmin_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.xmax = Label(segmentation_frame, text='x max size:', padx=2, pady=2)
        self.xmax.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.xmax_input = Entry(segmentation_frame, width=7)
        self.xmax_input.insert(0,'30')
        self.xmax_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.ymin = Label(segmentation_frame, text='y min size:', padx=2, pady=2)
        self.ymin.grid(row=3, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.ymin_input = Entry(segmentation_frame, width=7)
        self.ymin_input.insert(0,'2')
        self.ymin_input.grid(row=3, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.ymax = Label(segmentation_frame, text='y max size:', padx=2, pady=2)
        self.ymax.grid(row=4, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.ymax_input = Entry(segmentation_frame, width=7)
        self.ymax_input.insert(0,'30')
        self.ymax_input.grid(row=4, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.minMass = Label(segmentation_frame, text='mass min:', padx=2, pady=2)
        self.minMass.grid(row=5, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.minMass_input = Entry(segmentation_frame, width=7)
        self.minMass_input.insert(0,'0')
        self.minMass_input.grid(row=5, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.maxMass = Label(segmentation_frame, text='mass max:', padx=2, pady=2)
        self.maxMass.grid(row=6, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.maxMass_input = Entry(segmentation_frame, width=7)
        self.maxMass_input.insert(0,'10000')
        self.maxMass_input.grid(row=6, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.sigma = Label(segmentation_frame, text='blur sgima:', padx=2, pady=2)
        self.sigma.grid(row=7, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.sigma_input = Entry(segmentation_frame, width=7)
        self.sigma_input.insert(0,'0')
        self.sigma_input.grid(row=7, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.median = Label(segmentation_frame, text='median filter:', padx=2, pady=2)
        self.median.grid(row=8, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.median_input = Entry(segmentation_frame, width=7)
        self.median_input.insert(0,'0')
        self.median_input.grid(row=8, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.local = Label(segmentation_frame, text='local filter:', padx=2, pady=2)
        self.local.grid(row=9, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.local_input = Entry(segmentation_frame, width=7)
        self.local_input.insert(0,'0')
        self.local_input.grid(row=9, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        ROI_frame = LabelFrame(segmentation_frame, padx=2, pady=2, 
                                        width=100, text='ROI')
        ROI_frame.grid(row=11, column=0, columnspan=2, padx=5, pady=1, sticky='sew')
        
        self.ROI_x0y0 = Label(ROI_frame, text='x0, y0:', padx=2, pady=2)
        self.ROI_x0y0.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.ROI_x1y1 = Label(ROI_frame, text='x1, y1:', padx=2, pady=2)
        self.ROI_x1y1.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.ROIx0 = Entry(ROI_frame, width=7)
        self.ROIx0.insert(0,'0')
        self.ROIx0.grid(row=0, column=1)
        self.ROIy0 = Entry(ROI_frame, width=7)
        self.ROIy0.insert(0,'0')
        self.ROIy0.grid(row=0, column=2)
        self.ROIx1 = Entry(ROI_frame, width=7)
        self.ROIx1.insert(0,str(self.cam_res[0]))
        self.ROIx1.grid(row=1, column=1)
        self.ROIy1 = Entry(ROI_frame, width=7)
        self.ROIy1.insert(0,str(self.cam_res[1]))
        self.ROIy1.grid(row=1, column=2)
        
        
        
        segment_button = Button(segmentation_frame, text='Segment image', 
                                command = self.sementImage, padx=2, pady=7)
        segment_button.grid(row=12, column=0, padx=2, pady=7, sticky='ew')
        
        save_segment_button = Button(segmentation_frame, text='Save blobs', 
                                command = self.save_blobs, padx=2, pady=7)
        save_segment_button.grid(row=12, column=1, padx=2, pady=7, sticky='ew')
        
        
        
        
        
        
        
        # ====================================================================
        # (3) third column
        
        
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
        
        # save points button
        save_button = Button(button_frame, text='Save points', 
                                command = self.Save, padx=2, pady=4) 
        save_button.grid(row=6, column=0, padx=2, pady=2, sticky='ew')
        
        
        
        
        
        
        
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
        lab_frame.grid(row=2, column=0, columnspan=2, sticky='ew',
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
        
        
        
        
        # ================
        # sortgrid
        
        init_cal_frame = LabelFrame(Column3, text='4) sortgrid', 
                                  padx=2, pady=8, width=100)
        init_cal_frame.grid(row=1, column=0, columnspan=2, sticky='s', padx=2, 
                          pady=8)
        
        
        
        init_cal = Button(init_cal_frame, text='Initial calibration', 
                                command = self.calibrate, padx=2, pady=4) 
        init_cal.grid(row=0, column=0, padx=2, pady=2, sticky='new')
        
        
        show_cal = Button(init_cal_frame, text='Show calibration', 
                                command = self.show_calibration, padx=2, pady=4) 
        show_cal.grid(row=1, column=0, padx=2, pady=2, sticky='new')
        
        
        err_lf = Label(init_cal_frame)
        err_lf.grid(row=2, column=0)
        self.error = Label(err_lf, text='Error:', padx=2, pady=2)
        self.error.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.error_input = Label(err_lf, text='0.0', padx=2, pady=2,
                                 width=10, bg='white')
        self.error_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        sortgrid = Label(init_cal_frame)
        sortgrid.grid(row=3, column=0, pady=20)
        
        match_targets = Button(sortgrid, text='Match target file', 
                                command = self.matchTargetPoints, padx=2, pady=4) 
        match_targets.grid(row=0, column=0, padx=2, pady=2, sticky='new')
        
        
        save_cal_points = Button(sortgrid, text='Save cal points', 
                                command = self.saveTargetPoints, padx=2, pady=4) 
        save_cal_points.grid(row=1, column=0, padx=2, pady=2, sticky='new')
        
        
        
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
        self.cam = camera_Tsai(self.cam_name, cal_points_fname = cpf)
        self.cam.load('.')
        print('camera data loaded successfully.')
        cal = calibrate_Tsai(self.cam, self.cam.lab_points, self.cam.image_points)
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
        #image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
        #                     Image.ANTIALIAS)
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        
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
        self.cam.resolution = self.cam_res
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
        #image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
        #                     Image.ANTIALIAS)
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        
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
        #image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
        #                     Image.ANTIALIAS)
        
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        
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
        #image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
        #                     Image.ANTIALIAS)
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
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
    im_fname = '/home/ron/Desktop/Research/plankton_sweeming/experiments/20220916/MyPTV_analysis/Calibration/cal2.tif'    
    target_fname = '/home/ron/Desktop/Research/plankton_sweeming/experiments/20220916/MyPTV_analysis/Calibration/target_file'
    gui = initial_cal_gui('camX', im_fname, target_fname)

#cam = camera(cam_name, res)
#cam.save('.')

