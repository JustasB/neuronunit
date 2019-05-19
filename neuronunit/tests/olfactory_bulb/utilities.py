import numpy as np
import quantities as pq
import cPickle
import sys, os
import tempfile

class TestCache():
    pickle_dir = None
    cache = {}

    def __init__(self):
        self.pickle_dir = tempfile.gettempdir()

    def get(self, key):
        result = self.cache.get(key)

        if result is None:
            tmp_file =  os.path.join(self.pickle_dir, "TestCache." + str(key) + ".tmp")

            if os.path.exists(tmp_file):
                with open(tmp_file, "rb") as f:
                    result = cPickle.load(f)

            self.cache[key] = result

        return result

    def store(self, key, value):
        tmp_file = os.path.join(self.pickle_dir, "TestCache." + str(key) + ".tmp")

        with open(tmp_file, "wb") as f:
            cPickle.dump(value, f)

        self.cache[key] = value

    def clear(self):
        self.cache = {}

        tmp_files = [f for f in os.listdir(self.pickle_dir)
                     if f.startswith("TestCache.") and f.endswith(".tmp")]

        for f in tmp_files:
            os.remove(os.path.join(self.pickle_dir, f))


cache = TestCache()

def get_zero_crossings_neg2pos(voltage, after_delay=None):
    '''
    Returns the index locations where voltage value crossed 0 from neg->pos direction

    :param voltage: AnalogSignal or numpy array of voltage values
    :return: numpy array of 0-crossing indices
    '''
    if after_delay is not None:
        voltage = voltage.magnitude[np.where(voltage.times >= after_delay.rescale(pq.ms))]

    neg = voltage < 0
    return (neg[:-1] & ~neg[1:]).nonzero()[0]


def get_zero_crossings_pos2neg(voltage, after_delay=None):
    '''
    Returns the index locations where voltage value crossed 0 from pos->neg direction

    :param voltage: AnalogSignal or numpy array of voltage values
    :return: numpy array of 0-crossing indices
    '''

    if after_delay is not None:
        voltage = voltage.magnitude[np.where(voltage.times >= after_delay.rescale(pq.ms))]

    pos = voltage > 0
    return (pos[:-1] & ~pos[1:]).nonzero()[0]


def get_APs(voltage, ss_delay, method):
    '''

    :param voltage:
    :param ss_delay:
    :param method: 'd3v/dt3' or 'dv/dt=20'
    :return:
    '''
    # Find where voltage crosses 0 neg->pos
    crossings = get_zero_crossings_neg2pos(voltage)
    crossings = crossings[np.where(crossings > ss_delay * voltage.sampling_rate)]

    # Chop the series, keeping prior few ms
    pre_cross_window = voltage.sampling_rate * 5 * pq.ms
    cuts = np.array((crossings - pre_cross_window), dtype=int)
    ap_voltages = np.split(voltage, cuts)[1:]

    aps = []

    # Extract AP thresholds from the chopped series
    for i, v in enumerate(ap_voltages):
        if method == 'd3v/dt3':
            ap = extract_threshold_d3dt3(v)

        elif method == 'dv/dt=20':
            ap = extract_threshold_dvdt(v, 20)

        else:
            raise Exception("Unrecognized AP threshold extraction method " + str(method))

        # Only consider APs with detectable thresholds
        if ap is not None:
            aps.append(ap)

    # Chop again, but this time starting at the detected AP threshold
    cuts = np.array([ap["threshold_t"] * voltage.sampling_rate for ap in aps], dtype=int)
    ap_voltages = np.split(voltage, cuts)[1:]

    for i, ap in enumerate(aps):
        ap["voltage"] = ap_voltages[i]
        ap["peak_v"] = np.max(ap["voltage"])
        ap["amplitude"] = ap["peak_v"] - ap["threshold_v"]
        ap["half_width"] = len(np.where(ap["voltage"] > (ap["threshold_v"] + ap["amplitude"]/2.0))[0]) / ap["voltage"].sampling_rate

    return aps


def extract_threshold_d3dt3(v):
    # Compute the 3rd derivative
    v3 = np.diff(v.magnitude, n=3, axis=0)[:, 0]
    v3 = np.concatenate((v3, [0] * 3))

    # Zero-out very small fluctuations
    v3[np.where(np.abs(v3) < 1)] = 0

    # Get the first peak of v'''
    try:
        i_thresh = np.argmax(v3[0:get_zero_crossings_neg2pos(v3)[0] + 1])

        ap = {"threshold_v": v[i_thresh], "threshold_t": v.times[i_thresh].rescale(pq.ms)}

        return ap

    except:
        return None


def extract_threshold_dvdt(v, threshold):
    # Compute the 1st derivative
    v1 = np.diff(v.magnitude, n=1, axis=0)[:, 0]
    v1 = np.concatenate((v1, [0]))

    # Find the 20mV/ms crossing
    crossings = get_zero_crossings_neg2pos(v1 - threshold)

    if len(crossings) > 0:
        i_thresh = crossings[0]

        # Get the crossings value and time
        ap = {"threshold_v": v[i_thresh], "threshold_t": v.times[i_thresh].rescale(pq.ms)}

        return ap

    else:
        return None