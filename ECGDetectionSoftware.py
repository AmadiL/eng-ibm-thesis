#!/usr/bin/env python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.ttk import Frame, Button, Notebook, Treeview, Scrollbar
import scipy
from ECGAnalyzer import ECGAnalyzer
from matplotlib import style
from ttkthemes import ThemedStyle
import csv

SIGNALS = 'Sygnały'
SPECTRA = 'Widma'
INTERVALS = 'Interwały'
WAVES_TABLE = 'Tabela załamków'


class App(Frame):
    def __init__(self, master=None):
        super().__init__(master)
        style.use('bmh')  # matplotlib style
        self.width = 1280
        self.height = 720
        self.title = 'Wykrywanie charakterystycznych załamków w sygnale EKG - Amadeusz Lisiecki'
        self.ecg = None
        self.initUI()

    def initUI(self):
        self.style = ThemedStyle(self)
        self.style.set_theme("plastik")

        self.master.title(self.title)
        self.centerWindow()
        self.pack(fill=BOTH, expand=1)

        self.menubar = Menu(self.master)
        self.master.config(menu=self.menubar)

        self.fileMenu = Menu(self.menubar, tearoff=0)
        self.fileMenu.add_command(label="Otwórz próbkę", command=self.onOpenSample)
        self.fileMenu.add_command(label="Zapisz wszystko",  state=DISABLED, command=self.onSaveAll)
        self.fileMenu.add_command(label="Wyjdź", command=self.master.destroy)
        self.menubar.add_cascade(label="Plik", menu=self.fileMenu)

        self.topButtonFrame = Frame(self)
        self.topButtonFrame.pack(side=TOP, fill=X)
        self.openButton = Button(self.topButtonFrame, text="Otwórz próbkę", command=self.onOpenSample)
        self.openButton.pack(side=LEFT, padx=5, pady=5)
        self.saveAllButton = Button(self.topButtonFrame, state=DISABLED, text="Zapisz wszystko", command=self.onSaveAll)
        self.saveAllButton.pack(side=LEFT, padx=5, pady=5)
        self.file_label = Label(self.topButtonFrame, text='')
        self.file_label.pack(side=LEFT, padx=5, pady=5)

        self.menu_notebook = Notebook(self)
        self.menu_notebook.pack(fill=BOTH, expand=True)
        self.plots_notebook = Notebook(self.menu_notebook)
        self.plots_notebook.pack(fill=BOTH, expand=True)
        self.create_plot_tabs()
        self.menu_notebook.add(self.plots_notebook, text=SIGNALS)
        self.spectra_notebook = None
        self.table_frame = None
        self.table = None
        self.table_sb = None
        self.intervals_notebook = None

    def create_plot_tabs(self):
        # Plots
        if self.ecg is None:
            plots = ['Witaj!']
            self.plot_tabs = {}
        else:
            plots = ['Surowe EKG', 'Bez dryfu', 'Przefiltrowane', 'Zróżniczkowane', 'Do kwadratu', 'Scałkowane', 'Wynik']
            [tab['frame'].destroy() for tab in self.plot_tabs.values()]

        for i, tab in enumerate(plots):
            self.plot_tabs[tab] = {}
            self.plot_tabs[tab]['label'] = tab
            self.plot_tabs[tab]['frame'] = Frame(self.plots_notebook)
            self.plot_tabs[tab]['frame'].pack(fill=BOTH, expand=True)
            if tab != 'Witaj!':
                self.plot_tabs[tab]['btn'] = Button(self.plot_tabs[tab]['frame'], text='Zapisz dane', command=self.onSaveSample)
                self.plot_tabs[tab]['btn'].pack(side=TOP, padx=5, pady=5)
            self.plot_tabs[tab] = self.create_plot(self.plot_tabs[tab])
            self.draw_plot(self.plot_tabs[tab])
            self.plots_notebook.add(self.plot_tabs[tab]['frame'], text='{}. {}'.format(i + 1, tab))

    def create_spectra_tabs(self):
        # Spectra
        spectra = ['Surowe EKG', 'Przefiltrowane']
        if self.spectra_notebook is None:
            self.spectra_notebook = Notebook(self.menu_notebook)
            self.spectra_notebook.pack(fill=BOTH, expand=True)
            self.spectra_tabs = {}
            add_to_menu = True
        else:
            [tab['frame'].destroy() for tab in self.spectra_tabs.values()]
            add_to_menu = False

        for i, tab in enumerate(spectra):
            self.spectra_tabs[tab] = {}
            self.spectra_tabs[tab]['label'] = tab
            self.spectra_tabs[tab]['frame'] = Frame(self.spectra_notebook)
            self.spectra_tabs[tab]['frame'].pack(fill=BOTH, expand=True)
            self.spectra_tabs[tab] = self.create_plot(self.spectra_tabs[tab])
            self.draw_spectrum(self.spectra_tabs[tab])
            self.spectra_notebook.add(self.spectra_tabs[tab]['frame'], text='{}. {}'.format(i + 1, tab))

        if add_to_menu:
            self.menu_notebook.add(self.spectra_notebook, text=SPECTRA)

    def create_table(self):
        if self.table_frame is not None:
            self.table_frame.destroy()
        self.table_frame = Frame(self.menu_notebook)
        self.table_frame.pack(fill=BOTH, expand=True)
        self.table_save_button = Button(self.table_frame, text='Zapisz dane', command=self.onSaveSample)
        self.table_save_button.pack(side=TOP, padx=5, pady=5)
        self.table = Treeview(self.table_frame, show='headings', selectmode='browse')
        self.table['columns'] = tuple(sorted(self.ecg.indices.keys()) + ['Morfologia_T'])
        self.table.pack(side=BOTTOM, fill=BOTH, expand=True)

        self.table_sb = Scrollbar(self.table, orient="vertical", command=self.table.yview)
        self.table_sb.pack(side=RIGHT, fill='y')
        self.table.configure(yscrollcommand=self.table_sb.set)

        rows = 0
        for i, v in sorted(self.ecg.indices.items()) + [('Morfologia_T', self.ecg.T_shapes)]:
            rows = max(rows, len(v))
            text = i if i == 'Morfologia_T' else 'indeks {}'.format(i)
            self.table.heading(i, text=text)
            self.table.column(i, minwidth=0, width=100, stretch=True)
        for row in range(rows):
            text = []
            for i, v in sorted(self.ecg.indices.items()) + [('Morfologia_T', self.ecg.T_shapes)]:
                value = v[row] if row < len(v) else ''
                text.append(value)
            self.table.insert('', 'end', text=row, values=tuple(text))
        self.menu_notebook.add(self.table_frame, text=WAVES_TABLE)

    def create_intervals_tabs(self):
        # Intervals
        intervals = ['Wykres RR', 'Wykres QT', 'Tabela']
        PLOT_R = intervals[0]
        PLOT_QT = intervals[1]
        TABLE = intervals[2]
        if self.intervals_notebook is None:
            self.intervals_notebook = Notebook(self.menu_notebook)
            self.intervals_notebook.pack(fill=BOTH, expand=True)
            self.intervals_tabs = {}
            add_to_menu = True
        else:
            self.intervals_tabs[PLOT_R]['frame'].destroy()
            self.intervals_tabs[PLOT_QT]['frame'].destroy()
            self.intervals_tabs[TABLE].destroy()
            self.interval_table_sb.destroy()
            self.interval_table_frame.destroy()
            add_to_menu = False

        # Plot R
        self.intervals_tabs[PLOT_R] = {}
        self.intervals_tabs[PLOT_R]['label'] = PLOT_R
        self.intervals_tabs[PLOT_R]['frame'] = Frame(self.intervals_notebook)
        self.intervals_tabs[PLOT_R]['frame'].pack(fill=BOTH, expand=True)
        self.intervals_tabs[PLOT_R]['btn'] = Button(self.intervals_tabs[PLOT_R]['frame'], text='Zapisz dane',
                                                    command=self.onSaveSample)
        self.intervals_tabs[PLOT_R]['btn'].pack(side=TOP, padx=5, pady=5)
        self.intervals_tabs[PLOT_R] = self.create_plot(self.intervals_tabs[PLOT_R])
        self.draw_intervals(self.intervals_tabs[PLOT_R])
        self.intervals_notebook.add(self.intervals_tabs[PLOT_R]['frame'], text=PLOT_R)

        # Plot QT
        self.intervals_tabs[PLOT_QT] = {}
        self.intervals_tabs[PLOT_QT]['label'] = PLOT_QT
        self.intervals_tabs[PLOT_QT]['frame'] = Frame(self.intervals_notebook)
        self.intervals_tabs[PLOT_QT]['frame'].pack(fill=BOTH, expand=True)
        self.intervals_tabs[PLOT_QT]['btn'] = Button(self.intervals_tabs[PLOT_QT]['frame'], text='Zapisz dane',
                                                    command=self.onSaveSample)
        self.intervals_tabs[PLOT_QT]['btn'].pack(side=TOP, padx=5, pady=5)
        self.intervals_tabs[PLOT_QT] = self.create_plot(self.intervals_tabs[PLOT_QT])
        self.draw_intervals(self.intervals_tabs[PLOT_QT])
        self.intervals_notebook.add(self.intervals_tabs[PLOT_QT]['frame'], text=PLOT_QT)

        # Table
        self.interval_table_frame = Frame(self.menu_notebook)
        self.interval_table_frame.pack(fill=BOTH, expand=True)
        self.table_save_button = Button(self.interval_table_frame, text='Zapisz dane', command=self.onSaveSample)
        self.table_save_button.pack(side=TOP, padx=5, pady=5)
        self.intervals_tabs[TABLE] = Treeview(self.interval_table_frame, show='headings', selectmode='browse')
        self.intervals_tabs[TABLE]['columns'] = sorted(self.ecg.intervals.keys())
        self.intervals_tabs[TABLE].pack(side=BOTTOM, fill=BOTH, expand=True)

        self.interval_table_sb = Scrollbar(self.intervals_tabs[TABLE], orient="vertical", command=self.intervals_tabs[TABLE].yview)
        self.interval_table_sb.pack(side=RIGHT, fill='y')
        self.intervals_tabs[TABLE].configure(yscrollcommand=self.interval_table_sb.set)

        rows = 0
        for i, v in sorted(self.ecg.intervals.items()):
            rows = max(rows, len(v))
            text = '{} [ms]'.format(i) if i in ['QT', 'QTP', 'RR'] else i
            self.intervals_tabs[TABLE].heading(i, text=text)
            self.intervals_tabs[TABLE].column(i, minwidth=0, width=100, stretch=True)
        for row in range(rows):
            text = []
            for i, v in sorted(self.ecg.intervals.items()):
                value = v[row] if row < len(v) else ''
                text.append(value)
            self.intervals_tabs[TABLE].insert('', 'end', text=row, values=tuple(text))
        self.intervals_notebook.add(self.interval_table_frame, text=TABLE)

        if add_to_menu:
            self.menu_notebook.add(self.intervals_notebook, text=INTERVALS)

    def onOpenSample(self):
        try:
            ftypes = [('dat', '*.dat'), ('csv', '*.csv'), ('Wszystkie pliki', '*')]
            dlg = filedialog.Open(self, filetypes=ftypes)
            fl = dlg.show()

            if fl != '' and len(fl) != 0:
                self.file_label['text'] = 'Otwarty plik: {}'.format(fl)
                samplepath = re.sub(r'\.dat$', '', fl)
                self.ecg = ECGAnalyzer(samplepath)
                self.ecg.calculate()
                self.create_plot_tabs()
                self.create_spectra_tabs()
                self.create_intervals_tabs()
                self.create_table()
                self.menu_notebook.select(self.plots_notebook)
                self.saveAllButton['state'] = NORMAL
                self.fileMenu.entryconfig('Zapisz wszystko', state=NORMAL)
        except Exception as e:
            messagebox.showerror("Otwórz próbkę", e)

    def onSaveSample(self):
        try:
            dlg = filedialog.SaveAs(self)
            fl = dlg.show()
            if fl != '' and len(fl) != 0:
                if fl[-4:] != '.csv':
                    fl = '{}.csv'.format(fl)
                tab = self.menu_notebook.tab(self.menu_notebook.select(), 'text')
                if tab == SIGNALS:
                    signal_id = self.plots_notebook.index(self.plots_notebook.select())
                    self.save_signal(filename=fl, save_all=False, idn=signal_id)
                elif tab == SPECTRA:
                    pass
                elif tab == INTERVALS:
                    interval_tab_id = self.intervals_notebook.index(self.intervals_notebook.select())
                    self.save_interval(fl, interval_tab_id)
                elif tab == WAVES_TABLE:
                    self.save_waves_table(fl)
                else:
                    pass
        except Exception as e:
            messagebox.showerror("Zapisz dane", e)

    def onSaveAll(self):
        try:
            dlg = filedialog.SaveAs(self)
            fl = dlg.show()
            if fl != '' and len(fl) != 0:
                if fl[-4:] != '.csv':
                    fl = '{}.csv'.format(fl)
                self.save_all(fl)
        except Exception as e:
            messagebox.showerror("Zapisz wszystko", e)

    def save_waves_table(self, filename):
        indices_data = [[n, *i] for n, i in self.ecg.indices.items()]
        tm_data = ['T_morphology', *self.ecg.T_shapes]
        data = [*indices_data, tm_data]
        for d in data:
            append = max([len(i) for i in data]) - len(d)
            d += [''] * append
        data = zip(*data)
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            [csvwriter.writerow(row) for row in data]

    def save_interval(self, filename, idn):
        if idn == 0:
            data = ['RR', *self.ecg.intervals['RR']]
            data = [[i] for i in data]
        elif idn == 1:
            data = [[n, *i] for n, i in self.ecg.intervals.items() if n != 'RR']
            data = [*data]
        elif idn == 2:
            intervals_data = [[n, *i] for n, i in self.ecg.intervals.items()]
            data = [*intervals_data]
        if idn != 0:
            for d in data:
                append = max([len(i) for i in data]) - len(d)
                d += [''] * append
            data = zip(*data)
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            [csvwriter.writerow(row) for row in data]

    def save_signal(self, filename, save_all=False, idn=None):
        signals = [self.ecg.ecg_signal, self.ecg.no_drift_ecg_signal, self.ecg.filtered_ecg_signal,
                   self.ecg.differentiated_ecg_signal, self.ecg.squared_ecg_signal, self.ecg.integrated_ecg_signal,
                   self.ecg.filtered_ecg_signal]
        names = ['raw', 'no_drift', 'filtered', 'differentiated', 'squared', 'integrated', 'filtered']
        fs = ['fs', self.ecg.fs]
        if not save_all:
            name = ['title', names[idn]]
            signal = [names[idn], *signals[idn]]
            data = [name, fs, signal]
        else:
            name = ['title', 'all']
            signals_data = [[n, *s] for n, s in zip(names[:-1], signals[:-1])]
            indices_data = [[n, *i] for n, i in self.ecg.indices.items()]
            tm_data = ['T_morphology', *self.ecg.T_shapes]
            intervals_data = [[n, *i] for n, i in self.ecg.intervals.items()]
            data = [name, fs, *signals_data, *indices_data, tm_data, *intervals_data]

        for d in data:
            append = max([len(i) for i in data]) - len(d)
            d += [''] * append
        data = zip(*data)
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            [csvwriter.writerow(row) for row in data]

    def save_all(self, filename):
        self.save_signal(filename, save_all=True)

    def create_plot(self, tab):
        tab['fig'] = Figure(figsize=(3, 3), dpi=100)
        tab['plot'] = tab['fig'].add_subplot(111)
        tab['canvas'], tab['toolbar'] = self.create_plot_canvas(tab)
        return tab

    def create_plot_canvas(self, tab):
        canvas = FigureCanvasTkAgg(tab['fig'], tab['frame'])
        canvas.show()
        canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, tab['frame'])
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
        return canvas, toolbar

    def draw_intervals(self, tab):
        lines = []
        if tab['label'] == 'Wykres RR':
            title = 'Odstęp RR'
            x = self.ecg.t[self.ecg.indices['R'][1:]]
            y = self.ecg.intervals['RR']
            tab['plot'].plot(x, y)
            line = tab['plot'].plot(x, y, 'o')
            tab['plot'].set_title(title)
            tab['plot'].set_xlabel('Czas [s]')
            tab['plot'].set_ylabel('RR [ms]')
            tab['fig'].legend(line, ['RR'])
        else:
            title = 'Odstępy QT/QTP z i bez korekty'
            x = self.ecg.t[self.ecg.indices['T_end']]
            y = self.ecg.intervals['QT']
            lines.append(tab['plot'].plot(x, y, label='QT')[0])
            y = self.ecg.intervals['QTP']
            lines.append(tab['plot'].plot(x, y, label='QTP')[0])
            tab['plot'].set_xlabel('Czas [s]')
            tab['plot'].set_ylabel('QT/QTP [ms]')
            tab['plot'].tick_params('y')
            tab['plot2'] = tab['plot'].twinx()
            y = self.ecg.intervals['QTc']
            lines.append(tab['plot2'].plot(x, y, '--', label='QTc')[0])
            y = self.ecg.intervals['QTPc']
            lines.append(tab['plot2'].plot(x, y, '--', label='QTPc')[0])
            tab['plot2'].set_ylabel('skorygowany QT/QTP')
            tab['plot2'].tick_params('y')
            tab['fig'].tight_layout()
            tab['fig'].legend(lines, ['QT', 'QTP', 'QTc', 'QTPc'])
            tab['plot'].set_title(title)

    def draw_spectrum(self, tab):
        spectrum = None
        if tab['label'] == 'Surowe EKG':
            title = 'Widmo EKG'
            spectrum = self.ecg.spectrum(self.ecg.ecg_signal)
        elif tab['label'] == 'Przefiltrowane':
            title = 'Widmo EKG po filtracji'
            spectrum = self.ecg.spectrum(self.ecg.filtered_ecg_signal)
        if spectrum is not None:
            tab['plot'].set_title(title)
            line = tab['plot'].plot(spectrum.f, spectrum.y)
            tab['plot'].set_xlabel('Częstotliwość [Hz]')
            tab['plot'].set_ylabel('Amplituda')
            tab['fig'].legend(line, ['Widmo'])

    def draw_plot(self, tab):
        if self.ecg is None:
            self.draw_heart(tab['plot'])
        else:
            if tab['label'] != 'Wynik':
                x = self.ecg.t
                if tab['label'] == 'Surowe EKG':
                    title = 'EKG'
                    y = self.ecg.ecg_signal
                elif tab['label'] == 'Bez dryfu':
                    title = 'EKG + filtracja górnoprzepustowa'
                    y = self.ecg.no_drift_ecg_signal
                elif tab['label'] == 'Przefiltrowane':
                    title = 'EKG + filtracja górno- i dolnoprzepustowa'
                    y = self.ecg.filtered_ecg_signal
                elif tab['label'] == 'Zróżniczkowane':
                    title = 'EKG + filtracja + różniczkowanie'
                    y = self.ecg.differentiated_ecg_signal
                elif tab['label'] == 'Do kwadratu':
                    title = '(EKG + filtracja + różniczkowanie)²'
                    y = self.ecg.squared_ecg_signal
                elif tab['label'] == 'Scałkowane':
                    title = '(EKG + filtracja + różniczkowanie)² + całkowanie'
                    y = self.ecg.integrated_ecg_signal
                if y is not None and x is not None:
                    line = tab['plot'].plot(x, y, label=tab['label'])
                    tab['plot'].set_title(title)
                    tab['fig'].legend(line, [tab['label']])
            else:
                self.plot_result(tab)
            tab['plot'].axhline(0, color='red')
            tab['plot'].set_xlabel('Czas [s]')
            tab['plot'].set_ylabel('Amplituda')

    def plot_result(self, tab):
        title = 'EKG + filtracja + wynik detekcji załamków'
        lines = []
        tab['plot'].plot(self.ecg.t, self.ecg.filtered_ecg_signal)
        for i, v in sorted(self.ecg.indices.items()):
            x = self.ecg.t[v]
            y = self.ecg.filtered_ecg_signal[v]
            lines += tab['plot'].plot(x, y, 'o', label=i)
        tab['plot'].set_title(title)
        tab['fig'].legend(lines, sorted(self.ecg.indices.keys()))

    def draw_heart(self, a):
        x = scipy.linspace(-2, 2, 1000)
        y1 = scipy.sqrt(1 - (abs(x) - 1) ** 2)
        y2 = -3 * scipy.sqrt(1 - (abs(x) / 2) ** 0.5)
        a.fill_between(x, y1, color='red')
        a.fill_between(x, y2, color='red')
        a.set_xlim(-2.5, 2.5)
        a.text(0, -0.4, 'Analiza EKG:', fontsize=24, fontweight='bold',
               color='white', horizontalalignment='center')
        a.text(0, -0.8, 'Wykrywanie PQRST', fontsize=24, fontweight='bold',
               color='white', horizontalalignment='center')
        a.plot()

    def centerWindow(self):

        sw = self.master.winfo_screenwidth()
        sh = self.master.winfo_screenheight()

        x = (sw - self.width) / 2
        y = (sh - self.height) / 2
        self.master.geometry('%dx%d+%d+%d' % (self.width, self.height, x, y))


def main():
    try:
        root = Tk()
        app = App(master=root)
        app.mainloop()
    except Exception as e:
        messagebox.showerror('Error', e)


if __name__ == '__main__':
    main()
