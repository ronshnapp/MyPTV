# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A gui for manually marking points on a calibration target. 
"""


from PIL import Image, ImageTk
from tkinter import Label, Canvas, LabelFrame, Entry, Tk, Scrollbar, Button





class cal_point_gui(object):
    '''
    This is a Tkinter based graphical user interface that can be used to mark
    points on a static image, give their lab space coordinates, and then save
    the data as a text file in the format used by the calibration processes.
    '''
    
    
    def __init__(self, image_name, savename):
        '''
        input: 
        
        image_name - the calibration image used choose the points
        
        savename - the name of the text file in which the extracted data will
                   be saved.
        '''
        self.fname = savename
        self.image_name = image_name
        self.xy_marked = (-1, -1)
        self.point_list = []      # <-- list of points to save
        self.point_markers = []   # <-- position of crosses
        self.z = 1.0              # <-- Zoom level
        
        # set the window
        self.root = Tk()
        self.root.geometry('1000x700+320+70')
        self.root.title('MyPTV: calibration point marking GUI')
        
        # tracking_dic = {}
        # point_markers = []
        
        # load the image
        image = Image.open(self.image_name)
        image_size = image.size
        photo = ImageTk.PhotoImage(image)
        
        
        # buttons frame
        button_label = LabelFrame(self.root, padx=10, pady=10, width=100)
        button_label.grid(row=0, column=1, padx=(10), pady=10, sticky='nsew')
        
        # place the image inside a canvas in a frame
        self.board = Label(image=photo)
        frame = LabelFrame(self.root, bg='grey', width=700, height=450,
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
        
        # board.place(x=140, y=30)
        
        # setup click bottun
        self.canvas.bind("<Button-1>", self.location_handler)
        
        
        # selection indicators
        select_frame = LabelFrame(button_label, padx=2, pady=2)
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
        lab_frame = LabelFrame(button_label, padx=2, pady=2)
        lab_frame.grid(row=2, column=0, columnspan=2, sticky='ew',
                       padx=2, pady=2)
        self.x_input_label = Label(lab_frame, text='x lab:', padx=2, pady=2)
        self.y_input_label = Label(lab_frame, text='y lab:', padx=2, pady=2)
        self.z_input_label = Label(lab_frame, text='z lab:', padx=2, pady=2)
        self.x_input = Entry(lab_frame, width=14)
        self.x_input.insert(0,'0.0')
        self.y_input = Entry(lab_frame, width=14)
        self.y_input.insert(0,'0.0')
        self.z_input = Entry(lab_frame, width=14)
        self.z_input.insert(0,'0.0')
        
        self.x_input_label.grid(row=4, column=0, sticky='w', padx=2, pady=2)
        self.y_input_label.grid(row=5, column=0, sticky='w', padx=2, pady=2)
        self.z_input_label.grid(row=6, column=0, sticky='w', padx=2, pady=2)
        self.x_input.grid(row=4, column=1, sticky='w', padx=2, pady=2)
        self.y_input.grid(row=5, column=1, sticky='w', padx=2, pady=2)
        self.z_input.grid(row=6, column=1, sticky='w', padx=2, pady=2)
        

        # instructions tab
        instruction = LabelFrame(button_label, padx=2, pady=2,
                                 text='instructions')
        instruction.grid(row=0, columnspan=2, sticky='ew', padx=2, pady=2)
        text1 = 'Click on a calibration point, enter its lab coordinates in ' 
        text2 = 'the textboxes, and click Mark point. Repeat this for all the '
        text3 = 'relevant points, and click Save once finished to generate a ' 
        text4 = 'text file with points.\n\n'
        text5 = 'Zoom in/out with +/- \n\n'
        text6 = 'Shift-arrowkeys move the cursos in small increments.\n'
        text = text1 + text2 + text3 + text4 + text5 + text6
        instructions = Label(instruction, text=text, padx=2, pady=2, width=20,
                             wraplength=150, justify="left")
        instructions.grid(row=0,column=0)
        
        
        # set the mouse motion printing
        self.canvas.bind("<Motion>", self.motion)
        Mousepos_frame = LabelFrame(button_label, padx=2, pady=2)
        Mousepos_frame.grid(row=3, sticky='we', columnspan=2, padx=2, pady=2)
        self.Mousepos = Label(Mousepos_frame, text=' - , - ', padx=2, pady=2)
        self.Mousepos.grid(row=0, column=0, sticky='e', padx=2, pady=2)
        
        
        # zoom in and zoom out by pressing + and -
        self.root.bind('+', self.zoomIn)
        self.root.bind('-', self.zoomOut)
        
        # move cross with the arrow keys
        self.root.bind('<Shift-Left>', self.leftKey)
        self.root.bind('<Shift-Right>', self.rightKey)
        self.root.bind('<Shift-Up>', self.upKey)
        self.root.bind('<Shift-Down>', self.downKey)
        
        
        # Buttons frame
        button_frame = LabelFrame(button_label)
        button_frame.grid(row=4, column=0, columnspan=2, sticky='nwe', padx=2, 
                          pady=10)
        
        # add point button
        add_button = Button(button_frame, text='Mark point', 
                                command = self.addPoint, padx=2, pady=4)
        add_button.grid(row=0, column=0, padx=2, pady=2, sticky='ew')
        
        # forget last point button
        forget_last_button = Button(button_frame, text='Forget last point', 
                                command = self.forgetLast, padx=2, pady=4)
        forget_last_button.grid(row=1, column=0, padx=2, pady=2, sticky='ew')
        
        # save points button
        save_button = Button(button_frame, text='Save points', 
                                command = self.Save, padx=2, pady=4) 
        save_button.grid(row=2, column=0, padx=2, pady=2, sticky='ew')
        
        # quit button
        quit_button = Button(button_frame, text='Quit', 
                                command = self.Quit, padx=2, pady=4) 
        quit_button.grid(row=3, column=0, padx=2, pady=2, sticky='ew')
        
        
        
        # configure hte frames and run main loop
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        button_label.rowconfigure(0, weight=1)
        button_label.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        self.root.mainloop()
        
        
        
        


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
        self.z = self.z*1.1
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
    
        
    def zoomOut(self, event):
        '''zoom out with - key'''
        self.z = self.z*0.9
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
        '''save trajectory data'''
        from numpy import savetxt
        savetxt(self.fname, self.point_list, fmt='%.1f', delimiter='\t')
        print('Points saved at: %s'%self.fname)
        
        
    def Quit(self):
        '''quit the app'''
        self.root.destroy()




if __name__ == '__main__':
    im_fname = '/home/ron/Desktop/Research/plankton_sweeming/experiments/20220406/Cal/Cam2.tif'
    save_fname = 'test'
    gui = cal_point_gui(im_fname, save_fname)



