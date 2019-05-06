import numpy as np
import quantities as pq

def get_zero_crossings_neg2pos(voltage, after_delay=None):
    '''
    Returns the index locations where voltage value crossed 0 from neg->pos direction

    :param voltage: AnalogSignal or numpy array of voltage values
    :return: numpy array of 0-crossing indices
    '''
    if not after_delay:
        voltage = voltage.magnitude[np.where(voltage.times >= after_delay.rescale(pq.ms))]

    neg = voltage < 0
    return (neg[:-1] & ~neg[1:]).nonzero()[0]


def get_zero_crossings_pos2neg(voltage, after_delay=None):
    '''
    Returns the index locations where voltage value crossed 0 from pos->neg direction

    :param voltage: AnalogSignal or numpy array of voltage values
    :return: numpy array of 0-crossing indices
    '''

    if not after_delay:
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

    # Chop the series, keeping prior 10 ms
    pre_cross_window = voltage.sampling_rate * 10 * pq.ms
    cuts = np.array((crossings - pre_cross_window), dtype=int)
    ap_voltages = np.split(voltage, cuts)[1:]

    aps = []

    # Extract AP thresholds from the chopped series
    for i, v in enumerate(ap_voltages):
        if method == 'd3v/dt3':
            ap = extract_threshold_d3dt3(v)

        elif method == 'dv/dt=20':
            ap = extract_threshold_dvdt20(v)

        else:
            raise Exception("Unrecognized AP threshold extraction method " + str(method))

        aps.append(ap)

    cuts = np.array([ap["threshold_t"] * voltage.sampling_rate for ap in aps], dtype=int)
    ap_voltages = np.split(voltage, cuts)[1:]

    for i, ap in enumerate(aps):
        ap["voltage"] = ap_voltages[i]

    return aps


def extract_threshold_d3dt3(v):
    # Compute the 3rd derivative
    v3 = np.diff(v.magnitude, n=3, axis=0)[:, 0]
    v3 = np.concatenate((v3, [0] * 3))

    # Zero-out very small fluctuations
    v3[np.where(np.abs(v3) < 1)] = 0

    # Get the first peak of v'''
    i_thresh = np.argmax(v3[0:get_zero_crossings_neg2pos(v3)[0] + 1])

    ap = {"threshold_v": v[i_thresh], "threshold_t": v.times[i_thresh].rescale(pq.ms)}

    return ap


def extract_threshold_dvdt20(v):
    # Compute the 1st derivative
    v1 = np.diff(v.magnitude, n=1, axis=0)[:, 0]
    v1 = np.concatenate((v1, [0]))

    # Find the 20mV/ms crossing
    crossings = get_zero_crossings_neg2pos(v1 - 20)
    i_thresh = crossings[0]

    # Get the crossings value and time
    ap = {"threshold_v": v[i_thresh], "threshold_t": v.times[i_thresh].rescale(pq.ms)}

    return ap