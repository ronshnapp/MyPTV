# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A gui for the final calibration of a Tsai model camera. 
"""


#from PIL import Image, ImageTk
from tkinter import Label, LabelFrame, Entry, Tk, Button, Checkbutton, IntVar
from matplotlib.pyplot import subplots, show, imread


class cal_gui(object):
    '''
    This is a Tkinter based graphical user interface that can be used to mark
    points on a static image, give their lab space coordinates, and then save
    the data as a text file in the format used by the calibration processes.
    '''
    
    
    def __init__(self, calibrate_obj=None, cal_image=None):
        '''
        input: 
        
        calibrate_obj - An instance of the calibrate class, ready to perform
                        calibration.
        
        '''
        self.calibrate_obj = calibrate_obj
        self.cal_image = cal_image
        
        # set the window
        self.root = Tk()
        self.root.geometry('370x350')
        #self.root.resizable(0,0)
        self.root.resizable(height = None, width = None)
        self.root.title('MyPTV: Extended Zolof calibration GUI')


        # =============================
        # BUTTONS:
        
        # buttons frame
        button_label = LabelFrame(self.root, padx=10, pady=3, width=100)
        button_label.grid(row=0, column=0, padx=(10), pady=3, sticky='nsew')
        

        # Buttons frame1
        button_frame1 = Label(button_label)
        button_frame1.grid(row=0, column=0, columnspan=1, sticky='nsw', padx=2, 
                          pady=3)
        
        
        
        # external params calibration button
        calibrate_button = Button(button_frame1, text='Calibrate', 
                                command = self.calibrate, 
                                padx=2, pady=4, width=20)
        calibrate_button.grid(row=0, column=0, padx=10, pady=2, sticky='ew')
        
        
        # plot calibration
        plot_button = Button(button_frame1, text='Plot calibration', 
                                command = self.plot_calibration, 
                                padx=2, pady=4, width=20)
        plot_button.grid(row=1, column=0, padx=10, pady=2, sticky='ew')
        
        # plot error historam
        plot_button = Button(button_frame1, text='Plot error hist.', 
                                command = self.plot_err_hist, 
                                padx=2, pady=4, width=20)
        plot_button.grid(row=2, column=0, padx=10, pady=2, sticky='ew')
        
        
        # save points button
        save_button = Button(button_frame1, text='Save', 
                                command = self.Save, padx=10, pady=4, width=15) 
        save_button.grid(row=5, column=0, padx=10, pady=30, sticky='ew')
        
        
        # quit button
        quit_button = Button(button_frame1, text='Quit', 
                                command = self.Quit, padx=10, pady=4, width=15) 
        quit_button.grid(row=5, column=1, padx=10, pady=30, sticky='ew')
        
        
        
        
        
        
        
        # =============================
        # CALIBRATOR STATUS:
        
            
        status_label = LabelFrame(self.root, padx=10, pady=4, width=30)
        status_label.grid(row=1, column=0, padx=(10), pady=4, sticky='nsew')
        
        self.status = Label(status_label, text='Status:', padx=2, pady=2)
        self.status.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.status_show = Label(status_label, text='wating for action', padx=2, pady=2,
                                 width=30, anchor='w', fg='green')
        self.status_show.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        self.camera = Label(status_label, text='Camera:', padx=2, pady=2)
        self.camera.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        
        cam = 'camera name'
        if self.calibrate_obj is not None:
            cam = self.calibrate_obj.cam.name
        
        self.camera_show = Label(status_label, text=cam, padx=2, pady=2,
                                 width=30, anchor='w')
        self.camera_show.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)

        
        
        
        # =============================
        # CAMERA STATUS:
        
        
        
        dashboard = LabelFrame(self.root, padx=10, pady=10, width=100)
        dashboard.grid(row=2, column=0, padx=(10), pady=10, sticky='nsew')
        
        
        # err_Dashboard frame - where the error is shown
        err_dashboard = LabelFrame(dashboard)
        err_dashboard.grid(row=2, column=1, columnspan=1, sticky='sw', padx=2, 
                          pady=10)
        

        self.error = Label(err_dashboard, text='Error:', padx=2, pady=2)
        self.error.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.error_input = Label(err_dashboard, text='0.0', padx=2, pady=2,
                                 width=14, bg='white')
        self.error_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        err = self.calibrate_obj.mean_squared_err()
        self.error_input.config(text = '%.3e'%err)

        
        #if self.calibrate_obj is not None:
        #    self.update_cal_stats()
        
        
        
        
        # ==============================
        # RUN 
        
        # configure hte frames and run main loop
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        button_label.rowconfigure(0, weight=1)
        button_label.columnconfigure(0, weight=1)
        
        self.root.mainloop()
        
        
        
    def calibrate(self):
        '''Does the searchCalibration, i.e. external calibration'''
        print('\n','Iterating to minimize external parameters...','\n')
        self.status_show.configure(fg='red', 
                                   text='minimizing external parameters...')
        self.root.update()
        self.calibrate_obj.calibrate()
        err = self.calibrate_obj.mean_squared_err()
        print('\n','calibration error: %.3f pixels'%(err),'\n')
        self.error_input.config(text = '%.3e'%err)
        self.status_show.configure(fg='green', text='done! waiting for action')
    
    
        
    
    def plot_calibration(self):
        '''Plots the calibration using plot_proj function'''
        self.status_show.configure(fg='red', text='plotting calibration...')
        self.root.update()
        fig, ax = subplots()
        if self.cal_image is not None:
            img = imread(self.cal_image)
            ax.imshow(img, cmap='gray')
        self.calibrate_obj.plot_proj(ax=ax)
        show()
        self.status_show.configure(fg='green', text='done! waiting for action')
        
        
    def plot_err_hist(self):
        '''Plots the calibration error histogram'''
        self.status_show.configure(fg='red', text='plotting calibration...')
        self.root.update()
        fig, ax = subplots()
        self.calibrate_obj.plot_err_distribution(ax=ax)
        show()
        self.status_show.configure(fg='green', text='done! waiting for action')
        
     
    def Save(self):
        '''save the calibrated camera'''
        self.status_show.configure(fg='red', text='saving results...')
        self.root.update()
        self.calibrate_obj.cam.save('.')
        self.status_show.configure(fg='green', text='done! waiting for action')
        
        
    def Quit(self):
        '''quit the app'''
        self.status_show.configure(fg='red', text='bye!')
        self.root.update()
        self.root.destroy()




if __name__ == '__main__':
    gui = cal_gui()



