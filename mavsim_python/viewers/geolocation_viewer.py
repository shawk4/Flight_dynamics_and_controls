"""
target geolocation algorithm
    - Beard & McLain, PUP, 2012
    - Updated:
        4/1/2022 - RWB
        4/6/2022 - RWB
"""
import sys
sys.path.append('..')
from plotter.plotter import Plotter

class GeolocationViewer:
    def __init__(self, app,  dt = 0.01,
                 time_window_length = 30, # number of data points plotted at a time
                 plot_period = 0.2, # time interval between a plot update
                 data_recording_period = 0.1): # time interval between recording a data update
        self._dt = dt
        self._data_window_length= time_window_length/data_recording_period
        self._update_counter = 0
        self._plots_per_row = 4
        self._plotter = Plotter(app=app, plots_per_row=self._plots_per_row, 
                                window_width=1280, window_height=300)  # plot last time_window seconds of data
        self._plot_period = plot_period
        self._data_recording_period = data_recording_period
        self._plot_delay = 0
        self._data_recording_delay = 0
        self._time = 0
        red = (255,0,0)
        self._plotter.create_plot_widget(plot_id='err_n', xlabel='Time (s)', ylabel='err_n(m)',
                                        window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='err_e', xlabel='Time (s)', ylabel='err_e(m)',
                                        window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='err_d', xlabel='Time (s)', ylabel='err_d(m)',
                                        window_length=self._data_window_length)
        self._plotter.create_data_set(plot_id='err_n', data_label='err_n', data_color=red)
        self._plotter.create_data_set(plot_id='err_e', data_label='err_e', data_color=red)
        self._plotter.create_data_set(plot_id='err_d', data_label='err_d', data_color=red)
        self._plotter.show_window()

    def update(self, error):
        if self._data_recording_delay >= self._data_recording_period:
            self.__update_data(error, self._time)
            self._data_recording_delay = 0
        if self._plot_delay >= self._plot_period:
            self.__update_plot()
            self._plot_delay = 0
        self._plot_delay += self._dt
        self._data_recording_delay += self._dt
        self._time += self._dt

    def __update_data(self, error, t):
        self._plotter.add_data_point(plot_id='err_n', data_label='err_n', xvalue=t, yvalue=error.item(0))
        self._plotter.add_data_point(plot_id='err_e', data_label='err_e', xvalue=t, yvalue=error.item(1))
        self._plotter.add_data_point(plot_id='err_d', data_label='err_d', xvalue=t, yvalue=error.item(2))

    def __update_plot(self):
        self._plotter.update_plots()

    def close_geolocation_viewer(self):
        self._plotter.close_window()

    def save_plot_image(self, plot_name):
        self._plotter.save_image(plot_name)
