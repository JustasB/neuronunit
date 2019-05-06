import quantities

from neuronunit.tests.olfactory_bulb.utilities import get_zero_crossings_neg2pos

debug = True

from ..base import scores
import quantities as pq
from neuronunit.tests.olfactory_bulb import OlfactoryBulbCellTest, OlfactoryBulbCellSpikeTest
from matplotlib import pyplot as plt
from sciunit import capabilities as scap
from neuronunit import capabilities as ncap

from utilities import *
from abc import abstractmethod

class RheobaseTest(OlfactoryBulbCellTest):
    units = pq.nA
    score_type = scores.ZScore
    description = "A test of the rheobase current of a cell"
    name = "Rheobase current test"
    max_iterations = 20
    max_current = 10 * pq.nA
    min_precision = 0.005 # find the RB within this precision (e.g. 0.01 is 1%)

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.ProducesSpikes,)

    def generate_prediction(self, model):

        # if debug:
        #     import pydevd
        #     pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)

        model.set_stop_time(self.ss_delay + self.current_duration)

        # Get resting voltage
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": 0*pq.nA})

        resting_v = np.median(voltage.magnitude[-10:])

        # Get change in voltage in response to negative current
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": -0.1*pq.nA})

        response_v = np.median(voltage.magnitude[-10:])

        # Estimate the current necessary to generate APs from input resistance
        delta_v = response_v - resting_v
        nA_per_mV = -0.1 / delta_v
        nA_to_0mV = -resting_v * nA_per_mV * pq.nA

        upper_bound = nA_to_0mV
        lower_bound = 0 * pq.nA
        trial_current = lower_bound
        spikes_found = False
        iteration = 0

        if debug:
            trace_spike = None
            trace_nospike = None

        # Binary search to find the boundaries that enclose the true rheobase
        # Terminate when: spikes under no stim, or narrowed to 1%
        while iteration < self.max_iterations and \
              upper_bound > 0*pq.nA and \
              (upper_bound-lower_bound)/upper_bound > self.min_precision:


            model.set_stop_time(self.ss_delay + self.current_duration)

            voltage = model.inject_square_current({"delay":     self.ss_delay,
                                                   "duration":  self.current_duration,
                                                   "amplitude": trial_current})

            # Quickly check for APs without extracting other AP properties
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

            if len(crossings) >= 1:
                upper_bound = trial_current
                spikes_found = True

                if debug:
                    trace_spike = voltage
            else:
                lower_bound = trial_current

                if debug:
                    trace_nospike = voltage

            # Initially test if APs are produced within the maximum bounds
            if iteration > 0:
                trial_current = (upper_bound + lower_bound) / 2.0
            else:
                trial_current = upper_bound

            iteration += 1

        if debug:
            if trace_spike is not None:
                plt.plot(trace_spike.times, trace_spike)

            if trace_nospike is not None:
                plt.plot(trace_nospike.times, trace_nospike)

            plt.show()

        if upper_bound == 0 * pq.nA:
            raise Exception("Negative rheobase: Model produces spikes without stimulation")

        if not spikes_found:
            raise Exception("Rheobase was not found between 0 and " + str(self.nA_to_0mV))

        return upper_bound


class SagVoltageTest(OlfactoryBulbCellTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the sag voltage of a cell"
    name = "Sag voltage test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.SupportsVoltageClamp,)

    def generate_prediction(self, model):

        if debug:
            # import pydevd
            # pydevd.settrace('192.168.0.21', port=4200)
            pass

        model.set_temperature(self.temperature)

        # First, find the current that will bring the voltage to the testing voltage e.g. -90mV
        # Do this using a voltage clamp
        model.set_stop_time(self.ss_delay + self.current_duration)
        current = model.clamp_voltage(voltages= [self.sag_testing_voltage, self.sag_testing_voltage, 0] * pq.mV,
                                      durations=[self.ss_delay,            self.current_duration,    0]*pq.ms)

        inhibitory_current = np.median(current[-10:].magnitude)*pq.nA

        # Inject that current


        voltage = model.inject_square_current({"delay":self.ss_delay,
                                            "duration": self.current_duration,
                                           "amplitude": inhibitory_current})

        times = voltage.times

        # Find the minimum v within the ROI after SS and before end of injection
        t_roi = np.where((times > self.ss_delay) & (times <= self.ss_delay + self.current_duration))
        v_roi = voltage.base[t_roi]

        if debug:
            plt.plot(times[t_roi], v_roi)

        v_injection_ss = v_roi[-1]
        v_min_index = np.argmin(v_roi)
        v_min = v_roi[v_min_index]

        if debug:
            plt.plot([times[t_roi][-1]], [v_injection_ss], 'o', label="Current SS: "+str(v_injection_ss))
            plt.plot([times[t_roi][v_min_index]], [v_min], 'o', label="Min: "+str(v_min))
            plt.ylim((v_roi.min().magnitude, v_roi.max().magnitude))

        if v_min < v_injection_ss:
            sag = v_injection_ss - v_min

        # If there is no sag within the window, use the voltage at specific time point e.g. 100ms after current onset
        else:
            t_roi = np.where((times > self.ss_delay + self.sag_window - 1*pq.ms) & (times < self.ss_delay + self.sag_window + 1*pq.ms))
            v_roi = voltage.base[t_roi]
            v_at_sag_window = np.median(v_roi.magnitude) * pq.mV

            sag = v_injection_ss - v_at_sag_window

            if debug:
                plt.plot([np.median(times[t_roi])], [v_at_sag_window], 'o', label="Alt min: "+str(v_at_sag_window))


        if debug:
            plt.legend()
            plt.show()

        return sag


class ReboundSpikingTest(OlfactoryBulbCellSpikeTest):
    units = pq.dimensionless
    score_type = scores.BooleanScore
    description = "A test of the presence of rebound spikes"
    name = "Rebound spiking test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.SupportsVoltageClamp,
                             ncap.ProducesActionPotentials,)

    def generate_prediction(self, model):

        if debug:
            import pydevd
            pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)

        if self.rebound_ap_method == "sag":
            # First, find the current that will bring the voltage to the testing voltage e.g. -90mV
            # Do this using a voltage clamp
            model.set_stop_time(self.ss_delay + self.current_duration)
            current = model.clamp_voltage(voltages= [self.sag_testing_voltage, self.sag_testing_voltage, 0] * pq.mV,
                                          durations=[self.ss_delay,            self.current_duration,    0]*pq.ms)

            inhibitory_current = np.median(current[-10:].magnitude)*pq.nA

        elif self.rebound_ap_method == "-300pA":
            inhibitory_current = -0.3*pq.nA


        # Inject that current - but wait for additional time after end
        model.set_stop_time(self.ss_delay + self.current_duration + self.rebound_rest_time)

        voltage = model.inject_square_current({"delay":self.ss_delay,
                                                "duration": self.current_duration,
                                                "amplitude": inhibitory_current})

        # Look for APs after current is withdrawn
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay + self.current_duration)

        if debug:
            plt.plot(voltage.times, voltage)
            plt.show()

        return len(crossings) > 0


class AfterHyperpolarizationTest(OlfactoryBulbCellSpikeTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the sag voltage of a cell"
    name = "Sag voltage test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.ProducesSpikes,)

    def __init__(self, rheobase, *args, **kwargs):
        super(AfterHyperpolarizationTest, self).__init__(*args, **kwargs)
        self.rheobase = rheobase

    def get_first_ap(self, model):
        if debug:
            # import pydevd
            # pydevd.settrace('192.168.0.21', port=4200)
            pass

        model.set_temperature(self.temperature)

        # Inject rheobase
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": self.rheobase})

        aps = self.get_aps(voltage)

        return aps[0]

    def compute_amplitude(self, ap):
        if self.ahp_amplitude_method == 'threshold2min':
            i_min_v = np.argmin(ap["voltage"])

        elif self.ahp_amplitude_method == 'threshold2minWithin10ms':
            within10ms = np.where(ap["voltage"].times < ap["threshold_t"] + 10*pq.ms)
            i_min_v = np.argmin(ap["voltage"].magnitude[within10ms])

        else:
            raise Exception("Unrecognized AHP Amplitude method: " + str(self.ahp_amplitude_method))

        min_v = ap["voltage"][i_min_v]
        ahp_amplitude = ap["threshold_v"] - min_v

        return ahp_amplitude

    @abstractmethod
    def generate_prediction(self, model):
        pass


class AfterHyperpolarizationAmplitudeTest(AfterHyperpolarizationTest):
    def generate_prediction(self, model):

        ap = self.get_first_ap(model)

        if debug:
            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.show()

        return self.compute_amplitude(ap)


class AfterHyperpolarizationTimeTest(AfterHyperpolarizationTest):
    def generate_prediction(self, model):
        ap = self.get_first_ap(model)

        crossings = get_zero_crossings_pos2neg(ap["voltage"])
        crossing = crossings[0]
        crossing_time = ap["voltage"].times[crossing]

        if self.ahp_time_method == 'threshold2min':
            i_min_v = np.argmin(ap["voltage"])
            min_v_time = ap["voltage"].times[i_min_v]
            ahp_time = (min_v_time - crossing_time).rescale(pq.ms)

        elif self.ahp_time_method == 'threshold2amplitude50%':
            amplitude = self.compute_amplitude(ap)
            amp50 = ap["threshold_v"] - amplitude / 2.0
            i_amp50 = np.where(ap["voltage"].magnitude[:,0] < amp50)[0][0]
            t_amp50 = ap["voltage"].times[i_amp50]

            ahp_time = (t_amp50 - crossing_time).rescale(pq.ms)

        else:
            raise Exception("Unrecognized AHP Time method: " + str(self.ahp_time_method))

        if debug:
            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.show()

        return ahp_time


class FISlopeTest(OlfactoryBulbCellSpikeTest):
    units = pq.Hz/pq.nA
    score_type = scores.ZScore
    description = "A test of the FI curve slope/gain"
    name = "FI slope test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.ProducesActionPotentials,)

    def generate_prediction(self, model):

        if debug:
            import pydevd
            pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)
        model.set_stop_time(self.ss_delay + self.current_duration)

        currents = np.concatenate((np.arange(0, 300, 50), [300])) * pq.pA
        frequencies = []

        for i, trial_current in enumerate(currents):
            # Inject current
            voltage = model.inject_square_current({"delay":self.ss_delay,
                                                    "duration": self.current_duration,
                                                    "amplitude": trial_current})

            # Count APs
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

            freq = len(crossings) / self.current_duration.rescale(pq.sec)

            frequencies.append(freq)

        fi_gain = np.max(np.diff(frequencies))
        slope = fi_gain / 50 * 1000 * pq.Hz / pq.nA

        if debug:
            plt.plot(currents, frequencies, 'o')
            plt.show()

        return slope


class SpikeTrainTest(OlfactoryBulbCellSpikeTest):
    def __init__(self, rheobase, *args, **kwargs):
        super(SpikeTrainTest, self).__init__(*args, **kwargs)
        self.rheobase = rheobase

    def current_for_target_freq(self, model, rheobase, target_freq):
        model.set_stop_time(self.ss_delay + self.current_duration)

        freq_rb, _, _ = self.freq_at(model, rheobase)
        freq_2rb, _, _ = self.freq_at(model, rheobase * 2)
        freq_slope = (freq_2rb - freq_rb) / rheobase
        current_targetHz = rheobase + (target_freq - freq_rb) / freq_slope

        return current_targetHz

    def freq_at(self, model, current):
        # Inject current
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": current})

        # Count APs
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)
        freq = len(crossings) / self.current_duration.rescale(pq.sec)

        return freq, crossings, voltage


class ISICVTest(SpikeTrainTest):
    units = pq.dimensionless
    score_type = scores.ZScore
    description = "A test of the interspike interval (ISI) coefficient of variation (CV)"
    name = "ISI CV test"

    def generate_prediction(self, model):

        if debug:
            import pydevd
            pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)

        current_targetHz = self.current_for_target_freq(model,
                                                        self.rheobase,
                                                        self.spike_train_target_freq)

        freq_targetHz, crossings, voltage = self.freq_at(model, current_targetHz)

        isis = np.diff(crossings / voltage.sampling_rate)
        cv = np.std(isis) / np.mean(isis)

        if debug:
            plt.plot(isis, 'o')
            plt.show()

        return cv


class SpikeAccommodationTest(SpikeTrainTest):
    units = pq.dimensionless
    score_type = scores.ZScore
    description = "A test of the firing rate accommodation"
    name = "Spike Accommodation Test"

    def generate_prediction(self, model):

        if debug:
            import pydevd
            pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)

        if self.spike_train_method == "target_freq":
            current_targetHz = self.current_for_target_freq(model,
                                                            self.rheobase,
                                                            self.spike_train_target_freq)

            freq_targetHz, crossings, voltage = self.freq_at(model, current_targetHz)

        elif self.spike_train_method == "constant_current":
            model.set_stop_time(self.ss_delay + self.current_duration)

            voltage = model.inject_square_current({"delay": self.ss_delay,
                                                   "duration": self.current_duration,
                                                   "amplitude": self.spike_train_current})

            # Count APs
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)


        isis = np.diff(crossings / voltage.sampling_rate)
        isis.units = pq.sec

        ifr_first = 1.0/isis[0]
        ifr_last =  1.0/isis[-1]

        ifr_change = ifr_last - ifr_first
        ifr_change.units = pq.Hz

        if debug:
            plt.plot(isis.rescale(pq.ms), 'o')
            plt.show()

        return ifr_change


class SpikeAccommodationTimeConstantTest(SpikeTrainTest):
    units = pq.dimensionless
    score_type = scores.ZScore
    description = "A test of the firing rate accommodation time constant"
    name = "Spike Accommodation Time Constant Test"

    def generate_prediction(self, model):

        if debug:
            import pydevd
            pydevd.settrace('192.168.0.100', port=4200)

        model.set_temperature(self.temperature)

        if self.spike_train_method == "target_freq":
            current_targetHz = self.current_for_target_freq(model,
                                                            self.rheobase,
                                                            self.spike_train_target_freq)

            freq_targetHz, crossings, voltage = self.freq_at(model, current_targetHz)

        elif self.spike_train_method == "constant_current":
            model.set_stop_time(self.ss_delay + self.current_duration)

            voltage = model.inject_square_current({"delay": self.ss_delay,
                                                   "duration": self.current_duration,
                                                   "amplitude": self.spike_train_current})

            # Count APs
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        crossing_times = crossings / voltage.sampling_rate
        crossing_times -= crossing_times[0]

        isis = np.diff(crossing_times)
        isis.units = pq.sec

        ifrs = 1.0 / isis
        ifrs.units = pq.Hz
        ifrs = ifrs.magnitude
        ifr_times = (crossing_times - crossing_times[1])[1:]
        ifr_times = ifr_times.magnitude

        def ifr_func(t, start, finish, tau):
            return (start - finish) * np.exp(-t / tau) + finish

        from lmfit import Model

        model = Model(ifr_func)
        params = model.make_params(start=ifrs[0], finish=ifrs[-1], tau=10.0)
        params['tau'].min = 0
        result = model.fit(ifrs, t=ifr_times, params=params)

        start = result.best_values["start"]
        finish = result.best_values["finish"]
        tau = result.best_values["tau"]

        if debug:
            from matplotlib import pyplot as plt
            print(result.fit_report())

            plt.plot(crossing_times[1:], ifrs, 'bo')
            plt.plot(crossing_times[1:], result.best_fit, 'r-')
            plt.show()

        return tau * pq.ms



