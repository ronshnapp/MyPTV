# -*- coding: utf-8 -*-
"""
Created on April 23, 2022

@author: Ron

A gui for the final calibration of a Tsai model camera. 
"""


#from PIL import Image, ImageTk
from tkinter import (Label, LabelFrame, Entry, Tk, Button, Checkbutton,
                     IntVar, StringVar, OptionMenu, DISABLED, NORMAL)
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
        # NOTE: no hard-coded geometry() here on purpose. A fixed pixel
        # size clips the lower widgets whenever the actual required size
        # is larger than assumed - which happens with different fonts,
        # DPI/scaling, or Tk themes. Instead the window is sized from the
        # widgets' own requested size at the end of __init__, and that
        # size is enforced as a minimum so nothing can ever be cut off.
        self.root.resizable(True, True)
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
        # FORWARD-MODEL OPTIONS:
        #
        # [A] maps lab coordinates to pixels. Its three lab axes take
        # independent polynomial orders, because the depth direction is
        # usually far less distorted than the two in-plane ones. The
        # default (3,3,2) is the historical MyPTV basis.

        fwd_label = LabelFrame(self.root, text='Forward model [A] (lab -> pixels)',
                               padx=10, pady=4)
        fwd_label.grid(row=1, column=0, padx=(10), pady=4, sticky='nsew')

        self.a_order_x = StringVar(value='3')
        self.a_order_y = StringVar(value='3')
        self.a_order_z = StringVar(value='2')

        for col, (lbl, var) in enumerate([('X order:', self.a_order_x),
                                          ('Y order:', self.a_order_y),
                                          ('Z order:', self.a_order_z)]):
            Label(fwd_label, text=lbl, padx=2, pady=2).grid(
                row=0, column=2*col, sticky='nw', padx=2, pady=2)
            m = OptionMenu(fwd_label, var, '1', '2', '3')
            m.config(width=3)
            m.grid(row=0, column=2*col+1, sticky='nw', padx=2, pady=2)

        Label(fwd_label, anchor='w', justify='left', fg='gray25',
              text=('default (3,3,2) = the standard model.\n'
                    'the depth axis usually needs a lower order.')
              ).grid(row=1, column=0, columnspan=6, sticky='nw',
                     padx=2, pady=2)

        # =============================
        # BACKWARD-MODEL OPTIONS:
        #
        # These control how the camera-to-lab (epipolar line) model is
        # fit. The forward (lab-to-camera, [A]) model is unaffected by
        # them - which is why the "Error" readout below, which measures
        # only the forward projection, does NOT respond to these settings.
        # Watch "Ray err" instead when comparing them.

        opts_label = LabelFrame(self.root, text='Backward model (epipolar lines)',
                                padx=10, pady=4)
        opts_label.grid(row=2, column=0, padx=(10), pady=4, sticky='nsew')

        # -- polynomial order of B (ray directions) --
        # This applies to BOTH the fixed- and variable-origin models, so
        # it sits above the variable-origin checkbox and is never greyed.
        b_lbl = Label(opts_label, text='r(x) poly. order [B]:', padx=2, pady=2)
        b_lbl.grid(row=0, column=0, sticky='nw', padx=2, pady=2)

        self.b_order = StringVar(value='3')
        # order 4 is supported by calibrate(b_order=4) but deliberately
        # NOT offered here: with typical calibration point counts it
        # overfits badly, and the GUI should not hand that to a user who
        # has no reason to expect it.
        self.b_order_menu = OptionMenu(opts_label, self.b_order,
                                       '1', '2', '3')
        self.b_order_menu.config(width=12)
        self.b_order_menu.grid(row=0, column=1, sticky='nw', padx=2, pady=2)

        # -- variable origin on/off --
        self.var_origin = IntVar(value=0)
        self.var_origin_check = Checkbutton(
            opts_label, text='Use variable origin O(x)',
            variable=self.var_origin, command=self._toggle_origin_opts,
            padx=2, pady=2, anchor='w')
        self.var_origin_check.grid(row=1, column=0, columnspan=2,
                                   sticky='nw', padx=2, pady=2)

        # -- origin model: free / plane --
        self.origin_model_lbl = Label(opts_label, text='Origin model:',
                                      padx=2, pady=2)
        self.origin_model_lbl.grid(row=2, column=0, sticky='nw',
                                   padx=2, pady=2)

        self.origin_model = StringVar(value='plane')
        self.origin_model_menu = OptionMenu(opts_label, self.origin_model,
                                            'plane', 'free')
        self.origin_model_menu.config(width=12)
        self.origin_model_menu.grid(row=2, column=1, sticky='nw',
                                    padx=2, pady=2)

        # -- polynomial order of C --
        self.c_order_lbl = Label(opts_label, text='O(x) poly. order [C]:',
                                 padx=2, pady=2)
        self.c_order_lbl.grid(row=3, column=0, sticky='nw', padx=2, pady=2)

        self.c_order = StringVar(value='1')
        self.c_order_menu = OptionMenu(opts_label, self.c_order,
                                       '1', '2', '3')
        self.c_order_menu.config(width=12)
        self.c_order_menu.grid(row=3, column=1, sticky='nw', padx=2, pady=2)

        # -- hint --
        hint = Label(opts_label, anchor='w', justify='left', fg='gray25',
                     text=('lower order = less overfitting.\n'
                           '1=linear, 2=quadratic, 3=cubic\n'
                           '(3 for [B] reproduces the classic model)'))
        hint.grid(row=4, column=0, columnspan=2, sticky='nw', padx=2, pady=2)

        self._toggle_origin_opts()

        
        
        # =============================
        # CALIBRATOR STATUS:
        
            
        status_label = LabelFrame(self.root, padx=10, pady=4, width=30)
        status_label.grid(row=3, column=0, padx=(10), pady=4, sticky='nsew')
        
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
        dashboard.grid(row=4, column=0, padx=(10), pady=10, sticky='nsew')
        
        
        # err_Dashboard frame - where the error is shown
        err_dashboard = LabelFrame(dashboard)
        err_dashboard.grid(row=2, column=1, columnspan=1, sticky='sw', padx=2, 
                          pady=10)
        

        self.error = Label(err_dashboard, text='Error [px]:', padx=2, pady=2)
        self.error.grid(row=1, column=0, rowspan=1, sticky='nw', padx=2, pady=2)
        self.error_input = Label(err_dashboard, text='0.0', padx=2, pady=2,
                                 width=14, bg='white')
        self.error_input.grid(row=1, column=1, rowspan=1, sticky='nw', padx=2, pady=2)
        err = self.calibrate_obj.mean_squared_err()
        self.error_input.config(text = '%.3e'%err)

        # Ray (backward-model) error. Unlike "Error" above - which only
        # measures the forward [A] projection and is therefore identical
        # for every backward-model setting - this one actually responds to
        # the variable-origin options, so use it to compare them.
        self.ray_error = Label(err_dashboard, text='Ray err [lab]:',
                               padx=2, pady=2)
        self.ray_error.grid(row=2, column=0, rowspan=1, sticky='nw',
                            padx=2, pady=2)
        self.ray_error_input = Label(err_dashboard, text='--', padx=2, pady=2,
                                     width=14, bg='white')
        self.ray_error_input.grid(row=2, column=1, rowspan=1, sticky='nw',
                                  padx=2, pady=2)
        self._update_ray_err()

        
        #if self.calibrate_obj is not None:
        #    self.update_cal_stats()
        
        
        
        
        # ==============================
        # RUN 
        
        # configure the frames
        self.root.columnconfigure(0, weight=1)
        button_label.rowconfigure(0, weight=1)
        button_label.columnconfigure(0, weight=1)

        # Let every section keep its natural height, and give any EXTRA
        # vertical space to the dashboard at the bottom rather than
        # stretching the button block. (Previously row 0 alone had
        # weight=1, so the buttons absorbed all the slack.)
        for r in range(5):
            self.root.rowconfigure(r, weight=0)
        self.root.rowconfigure(4, weight=1)

        # Ask Tk to compute how much room the widgets actually need, then
        # open at exactly that size and refuse to shrink below it. This
        # keeps every button reachable no matter the font/DPI/theme.
        self.root.update_idletasks()
        req_w = self.root.winfo_reqwidth()
        req_h = self.root.winfo_reqheight()
        self.root.minsize(req_w, req_h)
        self.root.geometry('%dx%d' % (req_w, req_h))
        
        self.root.mainloop()



    def _toggle_origin_opts(self):
        '''
        Greys out the origin-model / polynomial-order selectors when the
        variable-origin model is switched off, since they have no effect
        in that case.
        '''
        state = NORMAL if self.var_origin.get() else DISABLED
        for w in [self.origin_model_lbl, self.origin_model_menu,
                  self.c_order_lbl, self.c_order_menu]:
            w.config(state=state)



    def _update_ray_err(self):
        '''
        Refreshes the ray-error readout. Before a calibration has been run
        the inlier lists don't exist yet, so this is a no-op then.
        '''
        try:
            rerr = self.calibrate_obj.mean_ray_err()
            self.ray_error_input.config(text='%.3e'%rerr)
        except Exception:
            self.ray_error_input.config(text='--')
        
        
        
    def calibrate(self):
        '''Does the searchCalibration, i.e. external calibration'''
        print('\n','Iterating to minimize external parameters...','\n')
        self.status_show.configure(fg='red', 
                                   text='minimizing external parameters...')
        self.root.update()

        use_var = bool(self.var_origin.get())
        kwargs = {'variable_origin': use_var,
                  'b_order': int(self.b_order.get()),
                  'a_order': (int(self.a_order_x.get()),
                              int(self.a_order_y.get()),
                              int(self.a_order_z.get()))}
        if use_var:
            kwargs['c_order'] = int(self.c_order.get())
            kwargs['origin_model'] = self.origin_model.get()

        print('calibrating with:', kwargs)

        # NOTE: deliberately NOT wrapped in a bare try/except that silently
        # falls back to a fixed-origin calibration. Such a fallback hides
        # real failures (too few surviving rays, bad correspondences, an
        # unsupported option) behind a normal-looking result, so you can
        # end up comparing settings that were never actually applied.
        # Failures are surfaced here instead.
        try:
            self.calibrate_obj.calibrate(**kwargs)
        except Exception as e:
            print('\nCALIBRATION FAILED: %s\n'%e)
            self.status_show.configure(fg='red',
                                       text='calibration FAILED - see console')
            return

        err = self.calibrate_obj.mean_squared_err()
        print('\n','calibration error: %.3f pixels'%(err),'\n')
        self.error_input.config(text = '%.3e'%err)
        self._update_ray_err()
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
