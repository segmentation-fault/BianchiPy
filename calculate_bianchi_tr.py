__author__ = 'antonio franco'

from scipy.optimize import fsolve
import matplotlib.pyplot as plt


class Bianchi:
    # Calculates the throughput of a saturated IEEE 802.11 WLAN basic scheme according to:
    # G.Bianchi, "Performance analysis of the IEEE 802.11 distributed coordination function," in IEEE
    # Journal on Selected Areas in Communications, vol. 18, no. 3, pp. 535 - 547, March 2000.
    # doi: 10.1109 / 49.840210
    def __init__(self, bitrate, n, ACK, SIFS, slot, DIFS, E_P, E_P_star, W, m, H, prop_delay):
        # INPUT:
        # bitrate: raw bitrate in bps
        # n: number of STAs
        # ACK: ACK length in bits
        # SIFS: SIFS duration in seconds
        # slot: slot duration in seconds
        # DIFS: DIFS duration in seconds
        # E_P: average packet payload size in bits
        # E_P_star:  average  length  of  the  longest  packetpayload involved in a collision in bit (for an example, see eq. 16 in the paper)
        # W: Minimum contention window size in slots
        # m: Retry limit
        # H: Header size in bits
        # prop_delay: Propagation delay in seconds
        #
        # OUTPUT:
        # S: normalized system throughput, defined as the fraction of time the channel is used to successfully transmit payload bits.
        # p: the probability of a collision seen by a packet being transmitted on the channel
        # t: the probability that a station transmits in a randomly chosen slot time
        # Ps: the probability that a transmission occurring on the channel is successful is given by the probability that exactly one station transmits on the channel, conditioned on the fact that at least one station transmits
        # P_tr: the probability that there is at least one transmission in the considered slot time
        # T_s: the average time the channel is sensed busy, in seconds
        # T_c: the average time the channel is sensed busy by each station during a collision in seconds.

        # independent varz
        self.bitrate = 56E6
        self.n = 30

        self.b00 = 0

        self.ACK = 14 * 8
        self.SIFS = 10E-6
        self.slot = 9E-6
        self.DIFS = 3 * 9E-6
        self.H = 2 * 8

        self.E_P = 80
        self.E_P_star = 100

        self.W = 12
        self.m = 7

        self.prop_delay = 0

        # dependent varz
        self.p = 0
        self.t = 0
        self.Ps = 0
        self.P_tr = 0
        self.T_s = 0
        self.T_c = 0
        self.S = 0

        self.bitrate = bitrate
        self.n = n
        self.ACK = ACK
        self.SIFS = SIFS
        self.slot = slot
        self.DIFS = DIFS
        self.E_P = E_P
        self.E_P_star = E_P_star
        self.W = W
        self.m = m
        self.H = H
        self.prop_delay = prop_delay

        self.calculate_p_t()
        self.calculate_Ptr()
        self.calculate_Ps()
        self.calculate_Ts()
        self.calculate_Tc()
        self.calculate_S()
        self.calculate_b00()

    def calculate_b00(self):
        self.b00 = 2.0 * (1.0 - 2.0 * self.p) * (1.0 - self.p) / ((1.0 - 2.0 * self.p) * (self.W + 1) + self.p * self.W * (1.0 - (2.0 * self.p)**self.m))

    def calculate_b(self, i, k):
        if i < self.m:
            W_i = 2.0 ** i * self.W
            b_i0 = self.p ** i * self.b00
        else:
            W_i = 2.0 ** self.m * self.W
            b_i0 = self.p ** self.m / (1 - self.p) * self.b00

        return (W_i - k) / W_i * b_i0

    def check_p_t(self):
        c1 = self.p - 1.0 + (1.0 - self.t) ** (self.n - 1.0) <= 1.49012e-08
        my_sum = 0.0
        for i in range(0, self.m):
            my_sum += (2.0 * self.p) ** i
        c2 = 2.0 / (1.0 + self.W + self.p * self.W * my_sum) - self.t <= 1.49012e-08
        return c1 and c2

    def calculate_p_t(self):
        def equations(x):
            p, t = x
            my_sum = 0.0
            for i in range(0, self.m):
                my_sum += (2.0 * p) ** i

            return ( p - 1.0 + (1.0 - t) ** (self.n - 1.0), 2.0 / (1.0 + self.W + p * self.W * my_sum) - t )

        self.p, self.t = fsolve(equations, (0.1, 0.1))

        print "System solved, error (p, tau): " + str(equations((self.p, self.t)))

    def calculate_Ps(self):
        self.Ps = self.n * self.t * (1.0 - self.t) ** (self.n - 1) / self.P_tr

    def calculate_Ptr(self):
        self.P_tr = 1.0 - (1.0 - self.t) ** self.n

    def calculate_Ts(self):
        self.T_s = self.H / self.bitrate + self.E_P / self.bitrate + self.SIFS + self.prop_delay + self.ACK / self.bitrate + self.DIFS + self.prop_delay

    def calculate_Tc(self):
        self.T_c = self.H / self.bitrate + self.E_P_star / self.bitrate + self.DIFS + self.prop_delay

    def calculate_S(self):
        self.S = self.Ps * self.P_tr * (self.E_P / self.bitrate) / (
        (1.0 - self.P_tr) * self.slot + self.P_tr * self.Ps * self.T_s + self.P_tr * (1.0 - self.Ps) * self.T_c)


class uniform_helper:
    # Calculates the parameters for a uniformly distributed packet length size between P_min and P_max (in bits)

    # independent varz
    P_min = 0
    P_max = 0

    #dependent varz
    E_P = 0
    E_P_star = 0

    def __init__(self, P_min, P_max):
        self.P_max = P_max
        self.P_min = P_min
        self.calculate_EP()
        self.calculate_E_P_star()

    def calculate_EP(self):
        self.E_P = (self.P_max - self.P_min) / 2.0

    def calculate_E_P_star(self):
        self.E_P_star = ( (self.P_max - self.P_min) ** 2.0 - 1 ) / (self.P_max - self.P_min)


if __name__ == "__main__":
    # Validation: see Fig.6 Bianchi paper

    bitrate = 1E6
    ACK = 112 + 128
    SIFS = 28E-6
    slot = 50E-6
    DIFS = 128E-6
    E_P = 8184
    E_P_star = E_P
    WW = [32, 128]
    mm = [3, 5]
    H = 272 + 128
    prop_delay = 0

    fig = plt.figure()

    ax = fig.gca()

    nn = range(5, 50)

    for W in WW:
        for m in mm:
            if m == 5 and W == 128:
                continue
            S = []
            for n in nn:
                B = Bianchi(bitrate, n, ACK, SIFS, slot, DIFS, E_P, E_P_star, W, m, H, prop_delay)
                S.append(B.S)

            if m == 3 and W == 128:
                marker = 'o'
            elif m == 3 and W == 32:
                marker = '^'
            else:
                marker = 's'

            plt.plot(nn, S, label="W = %d m = %d" % (W, m), marker=marker, c='k', markerfacecolor='w')

    ax.set_xlabel("Number of Stations")
    ax.set_ylabel("Saturation Throughput")

    leg = plt.legend(loc='best', fancybox=True, prop={'size': 12})
    frame = leg.get_frame()
    frame.set_alpha(0.5)  # make it semi-transparent

    plt.show()