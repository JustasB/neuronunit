import sciunit
import quantities as pq
from neo import AnalogSignal
import neuronunit
import numpy as np

class NeuronCellModel(sciunit.Model,
                      sciunit.capabilities.Runnable,
                      neuronunit.capabilities.ProducesMembranePotential,
                      neuronunit.capabilities.ReceivesSquareCurrent,
                      neuronunit.capabilities.SupportsVoltageClamp,
                      neuronunit.capabilities.ProducesSpikes,
                      neuronunit.capabilities.SupportsSettingTemperature):
    '''
    Defines a NeuronUnit model for running NeuronUnit tests against a
    cell model (1+ sections) implemented in NEURON simulator.

    The class implements methods to inject current and record membrane potential at specified cell segments.

    The class assumes the NEURON model has been loaded, synaptically isolated, and ready for current
    injection experiments. As input, it takes references to the NEURON segments where current is to be injected and
    membrane voltage measured.

    IMPORTANT: When modifying this class, ensure all unit tests pass before checking in your changes to prevent
    breaking of dependent NeuronUnit tests.

    Usage:

    # Load and setup your model in NEURON first
    from neuron import h
    h.load_file('cell.hoc')
    soma = h.Cell[0].soma
    dendrite = h.Cell[0].dend

    # Pass the segment where the current will be injected, and where the membrane potential will be measured
    model = NeuronCellModel(in_seg=soma(0.5), out_seg=dendrite(1.0), name="Smith et. al. (1996) Random Cell")

    # Judge the model
    test.judge(model)
    '''

    sampling_period = 0.25 # ms

    def __init__(self, in_seg, out_seg=None, name=None):
        super(NeuronCellModel, self).__init__()

        self.name = name
        self.in_seg = in_seg
        self.out_seg = out_seg if out_seg else in_seg

        from neuron import h
        self.h = h

        # Set up current and voltage clamps
        self.injector = self.h.IClamp(self.in_seg)
        self.vclamp = self.h.SEClamp(self.in_seg)
        self.reset_clamps()

        # Set up recorders for simulation time and membrane voltage
        self.tVector = self.h.Vector()
        self.vVector = self.h.Vector()
        self.vciVector = self.h.Vector()

        self.vVector.record(self.out_seg._ref_v, self.sampling_period)
        self.tVector.record(self.h._ref_t, self.sampling_period)
        self.vciVector.record(self.vclamp._ref_i, self.sampling_period)

    def get_backend(self):
        return self

    def set_stop_time(self, tstop):
        self.h.tstop = tstop.rescale(pq.ms).magnitude

    def set_temperature(self, celsius):
        self.h.celsius = celsius

    def inject_square_current(self, current = {"delay":0*pq.ms, "duration": 0*pq.ms, "amplitude": 0*pq.nA}):
        # Set the units that NEURON uses
        current["delay"].units = pq.ms
        current["duration"].units = pq.ms
        current["amplitude"].units = pq.nA

        self.injector.delay = float(current["delay"])
        self.injector.dur = float(current["duration"])
        self.injector.amp = float(current["amplitude"])

        self.h.run()
        self.reset_clamps()

        return self.get_membrane_potential()

    def clamp_voltage(self, voltages=[0*pq.mV, 0*pq.mV, 0*pq.mV], durations=[0]*3*pq.ms):
        self.vclamp.amp1, self.vclamp.amp2, self.vclamp.amp3 = [float(v) for v in voltages]
        self.vclamp.dur1, self.vclamp.dur2, self.vclamp.dur3 = [float(d) for d in durations]

        self.h.run()
        self.reset_clamps()

        return self.nrn_vector_to_AnalogSignal(self.vciVector, pq.nA)

    def get_membrane_potential(self):
        return self.nrn_vector_to_AnalogSignal(self.vVector, pq.mV)


    def nrn_vector_to_AnalogSignal(self, vector, units):
        '''
        Resample the signal stored by the NEURON vector at the specified steps_per_ms frequency

        :param vector: reference to a NEURON h.Vector()
        :param units: the units to use with the result
        :param steps_per_ms: the number of points to use to represent each ms of the recorded signal
        :return:
        '''
        #t = self.tVector.as_numpy()
        signal = vector.as_numpy()

        # new_t = np.linspace(t[0],t[-1],int(round(steps_per_ms*(t[-1]-t[0]))))
        # new_sig = np.interp(new_t, t, signal)

        return AnalogSignal(signal, sampling_period=self.sampling_period * pq.ms, units=units)


    def reset_clamps(self):
        '''
        Sets the current and voltage clamps to have no effect on the simulation
        '''

        self.injector.delay = 0
        self.injector.dur = 0
        self.injector.amp = 0

        # Set up voltage clamp
        self.vclamp.amp1, self.vclamp.amp2, self.vclamp.amp3 = [0,0,0]
        self.vclamp.dur1, self.vclamp.dur2, self.vclamp.dur3 = [0,0,0]
        self.vclamp.rs = 0.001  # MOhm