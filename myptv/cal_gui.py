# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A gui for manually marking points on a calibration target. 
"""


#from PIL import Image, ImageTk
from tkinter import Label, LabelFrame, Entry, Tk, Button
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
        self.err = self.calibrate_obj.mean_squared_err()
        
        # set the window
        self.root = Tk()
        self.root.geometry('400x380')
        self.root.resizable(0,0)
        self.root.title('MyPTV: calibration GUI')


        # =============================
        # BUTTONS:
        
        # buttons frame
        button_label = LabelFrame(self.root, padx=10, pady=3, width=100)
        button_label.grid(row=0, column=0, padx=(10), pady=3, sticky='nsew')
        

        # Buttons frame
        button_frame1 = Label(button_label)
        button_frame1.grid(row=0, column=0, columnspan=1, sticky='nsw', padx=2, 
                          pady=3)
        
        
        # external params calibration button
        external_button = Button(button_frame1, text='External calibration', 
                                command = self.external_calibration, 
                                padx=2, pady=4, width=20)
        external_button.grid(row=0, column=0, padx=10, pady=2, sticky='ew')
        
        # non-linear correction term
        fine_cal_button = Button(button_frame1, text='Fine calibration', 
                                command = self.fine_calibration, 
                                padx=2, pady=4, width=20)
        fine_cal_button.grid(row=1, column=0, padx=10, pady=2, sticky='ew')
        
        # plot calibration
        plot_button = Button(button_frame1, text='Plot calibration', 
                                command = self.plot_calibration, 
                                padx=2, pady=4, width=20)
        plot_button.grid(row=2, column=0, padx=10, pady=2, sticky='ew')
        
        
        # save points button
        save_button = Button(button_frame1, text='Save', 
                                command = self.Save, padx=10, pady=4, width=15) 
        save_button.grid(row=0, column=1, padx=10, pady=2, sticky='ew')
        
        # quit button
        quit_button = Button(button_frame1, text='Quit', 
                                command = self.Quit, padx=10, pady=4) 
        quit_button.grid(row=1, column=1, padx=10, pady=2, sticky='ew')
        
        
        
        
        
        
        
        # =============================
        # CALIBRATOR STATUS:
        
            
        status_label = LabelFrame(self.root, padx=10, pady=4, width=100)
        status_label.grid(row=1, column=0, padx=(10), pady=4, sticky='nsew')
        
        self.status = Label(status_label, text='Status:', padx=2, pady=2)
        self.status.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.status_show = Label(status_label, text='wating for action', padx=2, pady=2,
                                 width=50, anchor='w', fg='green')
        self.status_show.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        
        
        
        # =============================
        # CAMERA STATUS:
        
        
        
        dashboard = LabelFrame(self.root, padx=10, pady=10, width=100)
        dashboard.grid(row=2, column=0, padx=(10), pady=10, sticky='nsew')
        
        
        # pos_Dashboard frame - where the camera position is shown
        pos_dashboard = LabelFrame(dashboard)
        pos_dashboard.grid(row=1, column=0, columnspan=1, sticky='nsw', padx=10, 
                          pady=10)
        
        self.xloc = Label(pos_dashboard, text='Ox:', padx=2, pady=2)
        self.xloc.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.x_input = Entry(pos_dashboard, width=14)
        self.x_input.insert(0,'0.0')
        self.x_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.yloc = Label(pos_dashboard, text='Oy:', padx=2, pady=2)
        self.yloc.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.y_input = Entry(pos_dashboard, width=14)
        self.y_input.insert(0,'0.0')
        self.y_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.zloc = Label(pos_dashboard, text='Oz:', padx=2, pady=2)
        self.zloc.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.z_input = Entry(pos_dashboard, width=14)
        self.z_input.insert(0,'0.0')
        self.z_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        
        
        # ori_Dashboard frame - where the camera angles are shown
        ori_dashboard = LabelFrame(dashboard)
        ori_dashboard.grid(row=1, column=1, columnspan=1, sticky='nw', padx=2, 
                          pady=10)
        
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
        int_dashboard = LabelFrame(dashboard)
        int_dashboard.grid(row=2, column=0, columnspan=1, sticky='nw', padx=10, 
                          pady=10)
        
        self.f = Label(int_dashboard, text='f:', padx=2, pady=2)
        self.f.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.f_input = Entry(int_dashboard, width=14)
        self.f_input.insert(0,'0.0')
        self.f_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.xh = Label(int_dashboard, text='xh:', padx=2, pady=2)
        self.xh.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.xh_input = Entry(int_dashboard, width=14)
        self.xh_input.insert(0,'0.0')
        self.xh_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        self.yh = Label(int_dashboard, text='yh:', padx=2, pady=2)
        self.yh.grid(row=2, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.yh_input = Entry(int_dashboard, width=14)
        self.yh_input.insert(0,'0.0')
        self.yh_input.grid(row=2, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        
        # err_Dashboard frame - where the error is shown
        err_dashboard = LabelFrame(dashboard)
        err_dashboard.grid(row=2, column=1, columnspan=1, sticky='sw', padx=2, 
                          pady=10)
        
        self.error = Label(err_dashboard, text='Error:', padx=2, pady=2)
        self.error.grid(row=0, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.error_input = Label(err_dashboard, text='0.0', padx=2, pady=2,
                                 width=14, bg='white')
        self.error_input.grid(row=0, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        
        
        
        self.update_cal_stats()
        
        
        
        
        # ==============================
        # RUN 
        
        # configure hte frames and run main loop
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        button_label.rowconfigure(0, weight=1)
        button_label.columnconfigure(0, weight=1)
        
        self.root.mainloop()
        
        
        
        
        
    def update_cal_stats(self):
        ''' prints the calibrator status in the dashboard '''
        self.x_input.delete(0,'end')
        self.x_input.insert(0,'%.3f'%self.calibrate_obj.camera.O[0])
        
        self.y_input.delete(0,'end')
        self.y_input.insert(0,'%.3f'%self.calibrate_obj.camera.O[1])
        
        self.z_input.delete(0,'end')
        self.z_input.insert(0,'%.3f'%self.calibrate_obj.camera.O[2])
        
        self.xori_input.delete(0,'end')
        self.xori_input.insert(0,'%.4f'%self.calibrate_obj.camera.theta[0])
        
        self.yori_input.delete(0,'end')
        self.yori_input.insert(0,'%.4f'%self.calibrate_obj.camera.theta[1])
        
        self.zori_input.delete(0,'end')
        self.zori_input.insert(0,'%.4f'%self.calibrate_obj.camera.theta[2])
        
        self.f_input.delete(0,'end')
        self.f_input.insert(0,'%.1f'%self.calibrate_obj.camera.f)
        
        self.xh_input.delete(0,'end')
        self.xh_input.insert(0,'%.1f'%self.calibrate_obj.camera.xh)
        
        self.yh_input.delete(0,'end')
        self.yh_input.insert(0,'%.1f'%self.calibrate_obj.camera.yh)
        
        self.err = self.calibrate_obj.mean_squared_err()
        self.error_input.config(text = '%.3e px'%self.err)
        self.root.update()
        
        
    def external_calibration(self):
        '''Does the searchCalibration, i.e. external calibration'''
        print('\n','Iterating to minimize external parameters...','\n')
        self.status_show.configure(fg='red', 
                                   text='iterating to minimize external parameters...')
        self.root.update()
        self.calibrate_obj.searchCalibration(maxiter=2000)
        err = self.calibrate_obj.mean_squared_err()
        print('\n','calibration error: %.3f pixels'%(err),'\n')
        self.update_cal_stats()
        self.status_show.configure(fg='green', text='done! waiting for action')
    
    
    def fine_calibration(self):
        '''Does the fine (non-linear error) calibration'''
        print('\n', 'Iterating to minimize correction terms')
        self.status_show.configure(fg='red', 
                                   text='iterating to minimize fine calibration...')
        self.root.update()
        self.calibrate_obj.fineCalibration()
        err = self.calibrate_obj.mean_squared_err()
        print('\n','calibration error:', err,'\n')
        self.update_cal_stats()
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
        
        
    def Save(self):
        '''save the calibrated camera'''
        self.status_show.configure(fg='red', text='saving results...')
        self.root.update()
        self.calibrate_obj.camera.save('.')
        self.status_show.configure(fg='green', text='done! waiting for action')
        
        
    def Quit(self):
        '''quit the app'''
        self.status_show.configure(fg='red', text='bye!')
        self.root.update()
        self.root.destroy()




if __name__ == '__main__':
    gui = cal_gui()



