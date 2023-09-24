"""
TODO
add schematic picture
gain and freq. response calculations
"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline #interpolation
from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

DEFAULT_VALVE = 'E88CC'
DEFAULT_VSUPPLY = 265
DEFAULT_Ra = 33000
DEFAULT_Rk = 560
DEBUG = 1
DEFAULT_XMAX = 300
DEFAULT_YMAX = 10

class mclass:
    def __init__(self,  window):
        self.window = window
        self.fig = Figure(figsize=(13,8))
        self.ax = self.fig.add_subplot(111)

        # read valve specs and grid lines from csv
        self.read_valvedata()
        self.read_valvespecs()

        # build UI
        window.columnconfigure(0, weight=1)
        window.columnconfigure(1, weight=1)
        window.columnconfigure(2, weight=3)
        window.columnconfigure(3, weight=1)
        window.columnconfigure(4, weight=1)

        self.lbl_valve = Label(window, text="VALVE", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_valve.grid(row=1, column=0, sticky=W, padx=50, pady=5)
        self.str_valve = StringVar()
        self.str_valve.set(DEFAULT_VALVE)
        self.cmb_valve = ttk.Combobox(window, values=self.df['valve'].drop_duplicates().values.tolist(), textvariable=self.str_valve, font=('Courier New', 12), width=13)
        self.cmb_valve['state'] = 'readonly'
        self.cmb_valve.grid(row=1, column=1, sticky=W, padx=5, pady=5)
        self.cmb_valve.bind('<<ComboboxSelected>>', self.valve_changed)

        self.lbl_supply = Label(window, text="Vsupply, V", font=('Courier New', 12), background=self.window['bg'], width=20, anchor='w')
        self.lbl_supply.grid(row=2, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_supply = StringVar()
        self.str_supply.set(DEFAULT_VSUPPLY)
        self.etr_supply = Entry(window, textvariable=self.str_supply, font=('Courier New', 18), width=10)
        self.etr_supply.grid(row=2, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Ra = Label(window, text="Ra, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Ra.grid(row=3, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Ra = StringVar()
        self.str_Ra.set(DEFAULT_Ra)
        self.etr_Ra = Entry(window, textvariable=self.str_Ra, font=('Courier New', 18), width=10)
        self.etr_Ra.grid(row=3, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Rk = Label(window, text="Rk, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Rk.grid(row=4, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Rk = StringVar()
        self.str_Rk.set(DEFAULT_Rk)
        self.etr_Rk = Entry(window, textvariable=self.str_Rk, font=('Courier New', 18), width=10)
        self.etr_Rk.grid(row=4, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Vq = Label(window, text="Vq, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Vq.grid(row=5, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Vq = StringVar()
        self.lbl_Vqval = Label(window, textvariable=self.str_Vq, font=('Courier New', 18), width=10, anchor="w")
        self.lbl_Vqval.grid(row=5, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Iq = Label(window, text="Iq, mA", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Iq.grid(row=6, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Iq = StringVar()
        self.lbl_Iqval = Label(window, textvariable=self.str_Iq, font=('Courier New', 18), width=10, anchor="w")
        self.lbl_Iqval.grid(row=6, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Vgk = Label(window, text="Vgk, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Vgk.grid(row=7, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Vgk = StringVar()
        self.lbl_Vgkval = Label(window, textvariable=self.str_Vgk, font=('Courier New', 18), width=10, anchor="w")
        self.lbl_Vgkval.grid(row=7, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Xmax = Label(window, text="Xmax, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Xmax.grid(row=8, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Xmax = StringVar()
        self.str_Xmax.set(DEFAULT_XMAX)
        self.etr_Xmax = Entry(window, textvariable=self.str_Xmax, font=('Courier New', 18), width=10)
        self.etr_Xmax.grid(row=8, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Ymax = Label(window, text="Ymax, mA", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Ymax.grid(row=9, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Ymax = StringVar()
        self.str_Ymax.set(DEFAULT_YMAX)
        self.etr_Ymax = Entry(window, textvariable=self.str_Ymax, font=('Courier New', 18), width=10)
        self.etr_Ymax.grid(row=9, column=1, rowspan=1, sticky=W, padx=2, pady=5)


        # coordinates
        self.txt_coordinates = Text(bd=0, bg=window['bg'], fg='red', height=1, wrap="none", state="normal", font=('Courier New', 12), background=self.window['bg'])
        self.txt_coordinates.grid(row=12, column=1, columnspan=2, rowspan=1, sticky=W, padx=2, pady=5)
        self.txt_coordinates.config(highlightthickness = 0, borderwidth=0)
        self.txt_coordinates.config(state=DISABLED)

        # title
        self.title = Label(window, text='andmarti Loadline Plotter', fg='#1C5AAC', font=('Courier New', 24, 'bold'))
        self.title.grid(row=0, column=0, columnspan=5, padx=10, sticky=N)

        # plot
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().grid(row=1, column=4, rowspan=42, columnspan=3, sticky=W, padx=0, pady=40)
        self.canvas.mpl_connect('motion_notify_event', self.motion_hover)
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center', labelpad=20)
        self.ax.set_xlabel('Va, V', fontsize=16, loc='center',  labelpad=20)

        # bottom frame
        fm = Frame(window)
        self.lbl_cathodeloadline = Label(fm, text="Cathode loadline", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_cathodeloadline.grid(row = 0, column = 0, pady=0, sticky='w')
        self.chk_cathodeloadline_var = IntVar()
        self.chk_cathodeloadline = Checkbutton(fm, variable=self.chk_cathodeloadline_var, onvalue = 1, offvalue = 0, height=1, font=('Courier New', 12),
                                               command=self.chk_cathodeloadline_click, background=self.window['bg'], width=1, anchor="w")
        self.chk_cathodeloadline.select()
        self.chk_cathodeloadline.grid(row = 0, column = 1, columnspan=1, sticky='W')
        # valve specs
        self.str_specs = StringVar()
        self.str_specs.set("Specs: ")
        self.lbl_specs = Label(fm, textvariable=self.str_specs, font=('Courier New', 12), background=self.window['bg'])
        self.lbl_specs.grid(row = 1, column = 0, columnspan=62, sticky='W')

        #r = self.window.grid_size()[0]
        fm.grid(row=45, column=0, padx=2, pady=40, columnspan=10, sticky='W')

        #BUTTONS
        self.button_quit = Button(window, text="QUIT", command=self.quit, font=('Courier New', 18))
        self.button_quit.place(x=40, y=680)
        self.button_start = Button(window, text="PLOT", command=self.change_state, font=('Courier New', 18))
        self.button_start.place(x=160, y=680)
        self.button_clear = Button(window, text="CLEAR", command=self.clear_chart, font=('Courier New', 18), state='normal')
        self.button_clear.place(x=293, y=680)
        self.but_export = Button(window, text="EXPORT", command=self.export, font=('Courier New', 18))
        self.but_export.place(x=420, y=680)

        # bind focus out events
        self.etr_Xmax.bind("<FocusOut>", self.parameters_changed)
        self.etr_Ra.bind("<FocusOut>", self.parameters_changed)
        self.etr_supply.bind("<FocusOut>", self.parameters_changed)
        self.etr_Rk.bind("<FocusOut>", self.parameters_changed)
        self.etr_Ymax.bind("<FocusOut>", self.parameters_changed)

        self.etr_Xmax.bind("<Return>", self.parameters_changed)
        self.etr_Ra.bind("<Return>", self.parameters_changed)
        self.etr_supply.bind("<Return>", self.parameters_changed)
        self.etr_Rk.bind("<Return>", self.parameters_changed)
        self.etr_Ymax.bind("<Return>", self.parameters_changed)
        #end of ui
        # Read XMAX and YMAX from specs
        #self.updateMaxXY()
        self.valve_changed(None)

    def parameters_changed(self, event):
        self.change_state()

    def valve_changed(self, event):
        self.updateMaxXY()
        self.change_state()
        valve = self.str_valve.get()
        d=self.specs.loc[self.specs['valve']  == valve ]
        self.str_specs.set("Specs: mu: %s, ra: %s ohms, VaMax: %s V, Cga: %s pF, Cg to all but anode: %s pF, Pmax: %s W" % (
            d['mu'].iloc[0],
            d['ra'].iloc[0],
            d['VaMax'].iloc[0],
            d['Cga'].iloc[0],
            d['CgAEA'].iloc[0],
            d['Pmax'].iloc[0]
            ))

    # upodate xmax and ymax from valve specs
    def updateMaxXY(self):
        valve = self.str_valve.get()
        d=self.specs.loc[self.specs['valve']  == valve ]
        DEFAULT_XMAX = d['DefaultXmax'].iloc[0]
        DEFAULT_YMAX = d['DefaultYmax'].iloc[0]
        self.ax.set_ylim(0, d['DefaultYmax'].iloc[0])
        self.ax.set_xlim(0, d['DefaultXmax'].iloc[0])
        self.str_Ymax.set(DEFAULT_YMAX)
        self.str_Xmax.set(DEFAULT_XMAX)

    def clear_chart(self):
        self.ax.clear() # clear previous plot !!!!
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center', labelpad=20)
        self.ax.set_xlabel('Va, V', fontsize=16, loc='center',  labelpad=20)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def export(self):
        pass

    def chk_cathodeloadline_click(self):
        self.change_state()

    def quit(self):
        Tk().quit()

    def read_valvespecs(self):
        df = pd.read_csv('valves_specs.csv', comment='#', sep=',',  engine='python', skipinitialspace=True)
        #if DEBUG: print(df)
        self.specs = df

    def read_valvedata(self):
        pd.set_option('display.float_format', '{:20,.20f}'.format)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        #pd.set_option('display.width', None)
        #pd.set_option('display.max_colwidth', None)

        df = pd.read_csv('valves_data.csv', comment='#', sep='\r\n', names=['entire'], engine='python')
        df[['valve', 'curve', 'values']] = df['entire'].str.split(',', n=2, expand=True)
        df = df.drop(['entire'], axis=1)

        df['values'] = df['values'].apply(np.fromstring, sep=',')

        df = (df
          .assign(x=df['values'].str[::2], y=df['values'].str[1::2])
          .drop(columns='values')
          .explode(['x', 'y'])
        )
        df["x"] = pd.to_numeric(df["x"])
        df["y"] = pd.to_numeric(df["y"])

        df = df.set_index(['valve', 'curve' ])
        res = df.groupby(['valve', 'curve']).apply(lambda x: pd.Series(np.polyfit(x.x, x.y, 2), index=['a', 'b', 'c']))
        df = df.reset_index()
        self.df = pd.merge(df, res, how="inner", on=['valve', 'curve'])
        #if DEBUG: print(df)

    def motion_hover(self, event):
        if event.inaxes is not None:
            self.txt_coordinates.config(state=NORMAL)
            x = format(event.xdata, '.2f')
            y = format(event.ydata, '.2f')
            self.txt_coordinates.delete('1.0', END)
            self.txt_coordinates.insert(END, "(%s V, %s mA)" % (x, y))
            self.txt_coordinates.config(state=DISABLED)

    # function to plot loadline or refresh plot
    def change_state(self):
        # check values are valid
        if  self.can_convert_to_float(self.etr_Xmax.get()) == False: return
        if  self.can_convert_to_float(self.etr_Ra.get()) == False: return
        if  self.can_convert_to_float(self.etr_supply.get()) == False: return
        if  self.can_convert_to_float(self.etr_Rk.get()) == False: return
        if  self.can_convert_to_float(self.etr_Ymax.get()) == False: return

        self.clear_chart()

        # print(res)

        ymax = int(float(self.str_Ymax.get()))
        xmax = int(float(self.str_Xmax.get()))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        # select valve
        df = self.df
        valve = self.str_valve.get()
        de=df.loc[df['valve']  == valve ]

        # plot grid curves
        self.ax.set_title(valve, fontsize=24, pad=30, weight='bold')
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        #self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        #self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center')
        #self.ax.set_xlabel('Va, V', fontsize=16, loc='center')

        """
        #por regresion lineal
        num_pts = 100
        da=de[['valve','curve','a','b','c']]
        res = da.drop_duplicates()
        #print(res)
        for i in res.index:
            lower_limit = 0    #depende de curva
            x = np.linspace(lower_limit, xmax, num_pts)
            poly_coefs = res['a'][i], res['b'][i], res['c'][i]
            y = np.polyval(np.asarray(poly_coefs), x)
            #plt.plot(x, y, '-0')
            ind = np.argmin(y)
            x_min = x[ind]
            ax.set_ylim(0, ymax)
            ax.set_xlim(0, xmax)
            y_positive = y
            x_positive = x
            if res['curve'][i].strip() != 'Pmax':
                y_positive = y[x >= x_min]
                x_positive = x[x >= x_min]
            ax.plot(x_positive, y_positive, linewidth=2, label=res['curve'][i].strip())
        """

        # plot grid lines using existent points and interpolation
        self.ax.set_ylim(0, ymax)
        self.ax.set_xlim(0, xmax)
        da=de.loc[(de['curve'] != ' Pmax') & (de['x'] <= xmax) & (de['y'] <= ymax), ['valve','curve','x','y'] ]

        for name, g in da.groupby('curve'):
            #print(da)
            x_positive = g['x']
            y_positive = g['y']

            # try interpolation
            try:
                xnew = np.linspace(x_positive.min(), x_positive.max(), 300)
                spl = make_interp_spline(x_positive, y_positive, k=2)  # type: BSpline
                ynew = spl(xnew)
            except:
                # data points not enought ?
                xnew = x_positive
                ynew = y_positive
            self.ax.plot(xnew, ynew, '-', color='green', linewidth=2)

            # annotate grid curves
            if name == ' 0': name = "Vg = 0"
            self.ax.annotate(name + 'V', xy=(x_positive.max(), y_positive.max()), rotation=0, fontweight='bold')

        # plot Pmax using points
        da=de.loc[(de['curve'] == ' Pmax'), ['valve','curve','x','y'] ]
        x_positive = da['x']
        y_positive = da['y']
        # try interpolation
        xnew = np.linspace(x_positive.min(), x_positive.max(), 300)
        spl = make_interp_spline(x_positive, y_positive, k=2)  # type: BSpline
        ynew = spl(xnew)
        # add Pmax legend
        valve = self.str_valve.get()
        pmax=self.specs.loc[self.specs['valve']  == valve ]['Pmax'].iloc[0]
        self.ax.plot(xnew, ynew, 'r--', linewidth=2, label=da['curve'].iloc[0].strip() + ' = %sW' % str(pmax) )

        # plot loadline
        x_values = [float(self.str_supply.get()), 0]
        y_values = [0, float(self.str_supply.get()) * 1000 / float(self.str_Ra.get())]
        self.ax.plot(x_values, y_values, '-', color='blue', linewidth=2)

        # plot cathode line
        df = self.df
        valve = self.str_valve.get()
        de = df.loc[df['valve']  == valve ]
        de = de.drop(de[de['curve'] == ' Pmax'].index)
        data = de[['curve','a','b','c'] ].drop_duplicates()
        data[['curve', 'a', 'b', 'c']] = data[['curve', 'a', 'b', 'c']].apply(pd.to_numeric)
        data['y'] = - data['curve'] / float(self.str_Rk.get()) * 1000
        data['x'] = (-data['b']+(data['b']**2-4*data['a']*(data['c']-data['y']))**0.5)/(2*data['a'])
        x_positive = data['x']
        y_positive = data['y']
        # extrapolation
        xnew = np.linspace(x_positive.min(), x_positive.max(), 300)
        spl = make_interp_spline(x_positive, y_positive, k=1)  # type: BSpline
        ynew = spl(xnew)
        if self.chk_cathodeloadline_var.get() == 1:
            self.ax.plot(xnew, ynew, '-', color='orange', linewidth=2)
        # print(data)

        # quiscient
        b = float(self.str_supply.get()) * 1000 / float(self.str_Ra.get())
        a = -b / float(self.str_supply.get())
        yloadline = xnew * a + b
        idx = np.argwhere(np.diff(np.sign(ynew - yloadline))).flatten()
        self.vq = round(float(xnew[idx]), 2)
        self.str_Vq.set(format(self.vq, ".2f"))
        self.iq = round(float(yloadline[idx]), 2)
        self.str_Iq.set(format(self.iq, ".2f"))
        self.ax.plot(self.vq, self.iq, 'ro')

        # Vgk
        self.vgk = - self.iq / 1000 * float(self.str_Rk.get())
        self.str_Vgk.set(format(self.vgk, ".2f"))

        self.ax.legend(loc='upper left')
        self.canvas.draw()

    def can_convert_to_float(self, string):
        try:
            result = float(string)
            return True
        except ValueError:
            return False

window = Tk()
start = mclass(window)
window.mainloop()
