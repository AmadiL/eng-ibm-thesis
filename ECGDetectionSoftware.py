#!/usr/bin/env python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
from tkinter import filedialog
from tkinter.ttk import Frame, Button, Notebook, Treeview, Scrollbar
import scipy
from ECGAnalyzer import ECGAnalyzer
from matplotlib import style
from ttkthemes import ThemedStyle


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
        # self.style = Style()
        self.style = ThemedStyle(self)  # Like other Tkinter classes, a Style can take a master argument
        self.style.set_theme("plastik")
        # self.style.theme_use("arc")

        self.master.title(self.title)
        self.centerWindow()
        self.pack(fill=BOTH, expand=1)

        menubar = Menu(self.master)
        self.master.config(menu=menubar)

        fileMenu = Menu(menubar, tearoff=0)
        fileMenu.add_command(label="Otwórz próbkę", command=self.onOpenSample)
        menubar.add_cascade(label="Plik", menu=fileMenu)

        self.topButtonFrame = Frame(self)
        self.topButtonFrame.pack(side=TOP, fill=X)
        self.openButton = Button(self.topButtonFrame, text="Otwórz próbkę", command=self.onOpenSample)
        self.openButton.pack(side=LEFT, padx=5, pady=5)

        self.menu_notebook = Notebook(self)
        self.menu_notebook.pack(fill=BOTH, expand=True)
        self.plots_notebook = Notebook(self.menu_notebook)
        self.plots_notebook.pack(fill=BOTH, expand=True)
        self.create_plot_tabs()
        self.menu_notebook.add(self.plots_notebook, text='Sygnały')
        self.spectra_notebook = None
        self.table = None
        self.table_sb = None

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
            self.menu_notebook.add(self.spectra_notebook, text='Widma')

    def create_table(self):
        if self.table is not None:
            self.table.destroy()
            self.table_sb.destroy()
        self.table = Treeview(self.menu_notebook, show='headings', selectmode='browse')
        self.table['columns'] = tuple(sorted(self.ecg.indices.keys()))
        self.table.pack(side=LEFT)

        self.table_sb = Scrollbar(self.menu_notebook, orient="vertical", command=self.table.yview)
        self.table.pack(side=RIGHT, fill='y')

        self.table.configure(yscrollcommand=self.table_sb.set)

        rows = 0
        for i, v in sorted(self.ecg.indices.items()):
            rows = max(rows, len(v))
            self.table.heading(i, text=i)
            self.table.column(i, minwidth=0, width=100, stretch=True)
        for row in range(rows):
            text = []
            for i, v in sorted(self.ecg.indices.items()):
                value = v[row] if row < len(v) else None
                text.append(value)
            self.table.insert('', 'end', text=row, values=tuple(text))
        self.menu_notebook.add(self.table, text='Tabela')


    def onOpenSample(self):
        ftypes = [('dat', '*.dat'), ('Wszystkie pliki', '*')]
        dlg = filedialog.Open(self, filetypes=ftypes)
        fl = dlg.show()

        if fl != '' and len(fl) != 0:
            samplepath = re.sub(r'\.dat$', '', fl)
            self.ecg = ECGAnalyzer(samplepath)
            self.ecg.calculate()
            self.create_plot_tabs()
            self.create_spectra_tabs()
            self.create_table()

    def create_plot(self, tab):
        tab['fig'] = Figure(figsize=(5, 5), dpi=100)
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

    def draw_spectrum(self, tab):
        spectrum = None
        if tab['label'] == 'Surowe EKG':
            spectrum = self.ecg.spectrum(self.ecg.ecg_signal)
        elif tab['label'] == 'Przefiltrowane':
            spectrum = self.ecg.spectrum(self.ecg.filtered_ecg_signal)
        if spectrum is not None:
            tab['plot'].plot(spectrum.f, spectrum.y)
            tab['plot'].set_xlabel('Częstotliwość [Hz]')
            tab['plot'].set_ylabel('Amplituda')

    def draw_plot(self, tab):
        if self.ecg is None:
            self.draw_heart(tab['plot'])
        else:
            if tab['label'] != 'Wynik':
                x = self.ecg.t
                if tab['label'] == 'Surowe EKG':
                    y = self.ecg.ecg_signal
                elif tab['label'] == 'Bez dryfu':
                    y = self.ecg.no_drift_ecg_signal
                elif tab['label'] == 'Przefiltrowane':
                    y = self.ecg.filtered_ecg_signal
                elif tab['label'] == 'Zróżniczkowane':
                    y = self.ecg.differentiated_ecg_signal
                elif tab['label'] == 'Do kwadratu':
                    y = self.ecg.squared_ecg_signal
                elif tab['label'] == 'Scałkowane':
                    y = self.ecg.integrated_ecg_signal
                if y is not None and x is not None:
                    line = tab['plot'].plot(x, y, label=tab['label'])
                    tab['plot'].legend(line)
            else:
                self.plot_result(tab)
            tab['plot'].axhline(0, color='red')
            tab['plot'].set_xlabel('Czas [s]')
            tab['plot'].set_ylabel('Amplituda')
            tab['plot'].legend()

    def plot_result(self, tab):
        lines = []
        tab['plot'].plot(self.ecg.t, self.ecg.filtered_ecg_signal)
        tab['plot'].plot(self.ecg.t, self.ecg.differentiated_ecg_signal)
        for i, v in self.ecg.indices.items():
            x = self.ecg.t[v]
            y = self.ecg.filtered_ecg_signal[v]
            lines += tab['plot'].plot(x, y, 'o', label=i)
        tab['plot'].legend(lines)

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
    root = Tk()
    app = App(master=root)
    app.mainloop()


if __name__ == '__main__':
    main()
