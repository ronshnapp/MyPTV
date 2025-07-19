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

from myptv.imaging_mod import camera_wrapper, img_system
from myptv.segmentation_mod import particle_segmentation
from myptv.calibrate_mod import calibrate
from myptv.utils import match_calibration_blobs_and_points
from PIL import Image, ImageTk
from tkinter import Label,Canvas,LabelFrame,Entry,Tk,Scrollbar,Button,Listbox,END
from numpy import array, log10
import os



def fmt(a):
    if log10(abs(a))<3 and -3<log10(abs(a)):
        return '%.3f'%a
    else:
        return '%.3e'%a



class man_match_gui(object):
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
            cam = camera_wrapper(cn, './')
            cam.load()
            self.cam_list.append(cam)
            
        self.imsys = img_system(self.cam_list)
        
        
        # try to find a calibration folder
        ls = os.listdir('.')
        self.folder = '.'
        for fname in ls:
            if fname in ['calibration', 'Calibration', 'cal', 'Cal']:
                if os.path.isdir(os.path.join('.', fname)):
                    self.folder = os.path.join('.', fname)
        
        
        self.xy_marked = [(-1, -1) for i in range(len(self.cam_names))]
        self.point_list = []      # <-- list of points to save
        self.point_markers = []   # <-- position of crosses
        self.z = 1.0              # <-- Zoom level
        
        # set the window
        self.root = Tk()
        self.root.geometry('1100x700+320+70')
        self.root.title('MyPTV: Manual Matching GUI')
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
        # Second column
        
        
        Column3 = LabelFrame(self.root, padx=2, pady=10, width=100, 
                                  bg='#c2c2c2')
        Column3.grid(row=0, column=3, padx=(2), pady=10, sticky='nsew')
        
        
        #==================================
        # Marking points for initial calibration
        
        # Buttons frame
        button_frame = LabelFrame(Column3, text='Mark image points', 
                                  padx=2, pady=8, width=100)
        button_frame.grid(row=0, column=0, columnspan=2, sticky='nwe', padx=2, 
                          pady=8)
        
        
        
        select_frame = LabelFrame(button_frame, padx=2, pady=2)
        select_frame.grid(row=1, column=0, columnspan=2, sticky='ew',
                          padx=2, pady=2)
        
        forget_last_button = Button(select_frame, text='Forget point', 
                                command = self.forgetPoint, padx=2, pady=4)
        forget_last_button.grid(row=5, column=0, padx=2, pady=2, sticky='ew')
        
        
        
        
        # camera selector
        camera_select_frame = LabelFrame(button_frame, padx=2, pady=2)
        camera_select_frame.grid(row=0, column=0, columnspan=2, sticky='ew',
                          padx=2, pady=2)
        self.listbox = Listbox(camera_select_frame, selectbackground='blue')
        self.listbox.grid(row=0, column=0, rowspan=1, sticky='w', padx=2, pady=2)
        for cn in self.cam_names:
            self.listbox.insert(END, cn)
        
        
        self.listbox.bind("<<ListboxSelect>>", self.camera_select)
        self.listbox.select_set(0)
        self.current_camera_index = 0
        
        
        # selection indicators
        
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
                       padx=2, pady=30)
        
        calc_results_button = Button(lab_frame, text='Calculate in 3D', 
                                command = self.CalcRes, padx=2, pady=4)
        calc_results_button.grid(row=0, column=0, padx=2, pady=2, sticky='ew',
                                 columnspan=2)
        
        
        self.x_input_label = Label(lab_frame, text='x lab:', padx=2, pady=2)
        self.y_input_label = Label(lab_frame, text='y lab:', padx=2, pady=2)
        self.z_input_label = Label(lab_frame, text='z lab:', padx=2, pady=2)
        self.x_input = Entry(lab_frame, width=9)
        self.x_input.insert(0,'0.0')
        self.y_input = Entry(lab_frame, width=9)
        self.y_input.insert(0,'0.0')
        self.z_input = Entry(lab_frame, width=9)
        self.z_input.insert(0,'0.0')
        
        self.err_input_label = Label(lab_frame, text='error:', padx=2, pady=2)
        self.err_input = Entry(lab_frame, width=9)
        self.err_input.insert(0,'0.0')
        
        self.x_input_label.grid(row=4, column=0, sticky='w', padx=2, pady=2)
        self.y_input_label.grid(row=5, column=0, sticky='w', padx=2, pady=2)
        self.z_input_label.grid(row=6, column=0, sticky='w', padx=2, pady=2)
        self.err_input_label.grid(row=7, column=0, sticky='w', padx=2, pady=5)
        self.x_input.grid(row=4, column=1, sticky='w', padx=2, pady=2)
        self.y_input.grid(row=5, column=1, sticky='w', padx=2, pady=2)
        self.z_input.grid(row=6, column=1, sticky='w', padx=2, pady=2)
        self.err_input.grid(row=7, column=1, sticky='w', padx=2, pady=5)
        
        
        # set the mouse motion printing
        self.canvas.bind("<Motion>", self.motion)
        Mousepos_frame = Label(button_frame, padx=2, pady=2)
        Mousepos_frame.grid(row=9, sticky='we', columnspan=2, padx=2, pady=2)
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
        self.canvas.bind("<Button-4>", self.zoomIn)
        self.canvas.bind("<Button-5>", self.zoomOut)
        
        self.root.mainloop()

    

    
    def camera_select(self, event):
        selection = event.widget.curselection()
        self.current_camera_index = selection[0]
        print('Selected camera: ', self.cam_names[self.current_camera_index])
        
        image = self.images[self.current_camera_index]
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        
        photo = ImageTk.PhotoImage(image)
        
        self.board.configure(image = photo)
        self.board.image = photo
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=photo, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        x, y = self.xy_marked[self.current_camera_index]
        if x==-1 and y==-1:
            self.Xloc.configure(text = '-')
            self.Yloc.configure(text = '-')
        else:
            self.Xloc.configure(text = x)
            self.Yloc.configure(text = y)
        
        self.mark_points()
    


    def CalcRes(self):
        '''add marked point to list'''
        coords = {}
        for i in range(len(self.cam_names)):
            if self.xy_marked[i][0] != -1 and self.xy_marked[i][1] != -1:
                coords[i] = self.xy_marked[i]
        res = self.imsys.stereo_match(coords, d_max = 1e10)
        print(res)
        
        xyz = res[0]
        err = res[-1]
        print(( fmt(xyz[0]), fmt(xyz[1]), fmt(xyz[2])))
        
        self.x_input.delete(0, END)
        self.x_input.insert(0, fmt(xyz[0]))
        self.y_input.delete(0, END)
        self.y_input.insert(0, fmt(xyz[1]))
        self.z_input.delete(0, END)
        self.z_input.insert(0, fmt(xyz[2]))
        self.err_input.delete(0, END)
        self.err_input.insert(0, fmt(err))
            
        
    def forgetPoint(self):
        self.xy_marked[self.current_camera_index] = (-1, -1)
        self.Xloc.configure(text = '-')
        self.Yloc.configure(text = '-')
        self.mark_points()
        
        
    def rightKey(self, event):
        '''right key = move cross'''
        i = self.current_camera_index
        self.xy_marked[i] = (self.xy_marked[i][0]+1, self.xy_marked[i][1])
        self.mark_points()
        self.Xloc.configure(text = self.xy_marked[i][0]) 
        
        
    def leftKey(self, event):
        '''left key = move cross'''
        i = self.current_camera_index
        self.xy_marked[i] = (self.xy_marked[i][0]-1, self.xy_marked[i][1])
        self.mark_points()
        self.Xloc.configure(text = self.xy_marked[i][0])
        
    
    def downKey(self, event):
        '''left key = move cross'''
        i = self.current_camera_index
        self.xy_marked[i] = (self.xy_marked[i][0], self.xy_marked[i][1]+1)
        self.mark_points()
        self.Yloc.configure(text = self.xy_marked[i][1])
        
        
    def upKey(self, event):
        '''left key = move cross'''
        i = self.current_camera_index
        self.xy_marked[i] = (self.xy_marked[i][0], self.xy_marked[i][1]-1)
        self.mark_points()
        self.Yloc.configure(text = self.xy_marked[i][1])
        
    
    
    def zoomIn(self, event):
        '''zoom in the image with + key'''
        self.z = self.z*1.15
        image = self.images[self.current_camera_index]
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        self.mark_points()
        
    
        
    def zoomOut(self, event):
        '''zoom out with - key'''
        self.z = self.z*(1/1.15)
        image = self.images[self.current_camera_index]
        s = image.size
        image = image.resize((int(s[0]*self.z),int(s[1]*self.z)),
                             Image.LANCZOS)
        new_bird = ImageTk.PhotoImage(image)
        self.board.configure(image = new_bird)
        self.board.image = new_bird
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, image=new_bird, anchor='nw')
        self.canvas.configure(scrollregion = self.canvas.bbox("all"))
        
        self.mark_points()
        
        
        
    def location_handler(self, event):
        '''handling the location of the mouse pointer'''
        #x,y = int(event.x/z), int(event.y/z)
        
        x = int(self.canvas.canvasx(event.x)/self.z) + int(self.hbar.get()[1])
        y = int(self.canvas.canvasy(event.y)/self.z) + int(self.vbar.get()[1])
        
        self.Xloc.configure(text = x) 
        self.Yloc.configure(text = y)
        self.xy_marked[self.current_camera_index] = (x, y)
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

        x, y = self.xy_marked[self.current_camera_index]
        x_ = int(x*self.z) - int(self.hbar.get()[1])
        y_ = int(y*self.z) - int(self.vbar.get()[1])
        c1 = self.canvas.create_line(x_-5, y_, x_+5, y_, 
                                     fill="#e31010", width=2)
        c2 = self.canvas.create_line(x_, y_-5, x_, y_+5, 
                                     fill="#e31010", width=2)    
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


    def Quit(self):
        '''quit the app'''
        self.root.destroy()




if __name__ == '__main__':
    im_fname = ['../example/Calibration/cal1.tif',
                '../example/Calibration/cal2.tif',
                '../example/Calibration/cal3.tif']
    
    cameras_folder = '../example'
    
    camera_names = ['cam1', 'cam2', 'cam3']
    
    gui = man_match_gui(camera_names, im_fname, cameras_folder=cameras_folder)


