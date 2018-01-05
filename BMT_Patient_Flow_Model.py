import simpy
from numpy.random import RandomState

"""
Simple BMT Patient Flow Model
based on: © 2017 Mark Isken
http://hselab.org/simpy-first-oo-patflow-model.html
https://github.com/misken/obsimpy

"""

# Arrival rate and length of stay inputs.
ARR_RATE = 0.2
MEAN_LOS_RFL = .03
MEAN_LOS_CN1 = .05
MEAN_LOS_PTP = .5
MEAN_LOS_CN2 = .05
MEAN_LOS_CTH = .5
MEAN_LOS_STM = .5
MEAN_LOS_TPT = .8
MEAN_LOS_AFO = .05
MEAN_LOS_MTC = .3
MEAN_LOS_AFG = .25
MEAN_LOS_LTF = .005

# Unit capacities
CAPACITY_RFL = 2
CAPACITY_CN1 = 6
CAPACITY_PTP = 24
CAPACITY_CN2 = 5
CAPACITY_CTH = 10
CAPACITY_STM = 5
CAPACITY_TPT = 5
CAPACITY_AFO = 10
CAPACITY_MTC = 10
CAPACITY_AFG = 25
CAPACITY_LTF = 10


# Random number seed
RNG_SEED = 3590


class BMT_unit(object):
    """ Models a BMT_unit with fixed capacity.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        name : str
            unit name
        capacity : integer (or None)
            Number of beds. Use None for infinite capacity.

    """

    def __init__(self, env, name, capacity=None, debug=False):
        if capacity is None:
            self.capacity = simpy.core.Infinity
        else:
            self.capacity = capacity

        self.env = env
        self.name = name
        self.debug = debug

        # Use a simpy Resource as one of the class members
        self.unit = simpy.Resource(env, capacity)

        # Statistical accumulators
        self.num_entries = 0
        self.num_exits = 0
        self.tot_occ_time = 0.0

        # The out member will get set to destination unit
        self.out = None

    def put(self, bmtp):
        """ A process method called when a bed is requested in the unit.

            The logic of this method is reminiscent of the routing logic
            in the process oriented obflow models 1-3. However, this method
            is used for all of the units - no repeated logic.

            Parameters
            ----------
            env : simpy.Environment
                the simulation environment

            bmtp : BMT Patient object
                the patient requesting the bed

        """

        if self.debug:
            print("{} trying to get {} at {:.4f}".format(bmtp.name,
                                                    self.name, self.env.now))

        # Increments patient's attribute number of units visited
        bmtp.current_stay_num += 1
        # Timestamp of request time
        bed_request_ts = self.env.now
        # Request a bed
        bed_request = self.unit.request()
        # Store bed request and timestamp in patient's request lists
        bmtp.bed_requests[bmtp.current_stay_num] = bed_request
        bmtp.request_entry_ts[bmtp.current_stay_num] = self.env.now
        # Yield until we get a bed
        yield bed_request

        # Seized a bed.
        bmtp.entry_ts[bmtp.current_stay_num] = self.env.now

        # Check if we have a bed from a previous stay and release it.
        # Update stats for previous unit.

        if bmtp.bed_requests[bmtp.current_stay_num - 1] is not None:
            previous_request = bmtp.bed_requests[bmtp.current_stay_num - 1]
            previous_unit_name = \
                bmtp.planned_route_stop[bmtp.current_stay_num - 1]
            previous_unit = bmtp_units[previous_unit_name]
            previous_unit.unit.release(previous_request)
            previous_unit.num_exits += 1
            previous_unit.tot_occ_time += \
                self.env.now - bmtp.entry_ts[bmtp.current_stay_num - 1]
            bmtp.exit_ts[bmtp.current_stay_num - 1] = self.env.now

        if self.debug:
            print("{} entering {} at {:.4f}".format(bmtp.name, self.name,
                                                self.env.now))
        self.num_entries += 1
        if self.debug:
            if self.env.now > bed_request_ts:
                waittime = self.env.now - bed_request_ts
                print("{} waited {:.4f} time units for {} bed".format(bmtp.name,
                                                                  waittime,
                                                                  self.name))

        # Determine los and then yield for the stay
        los = bmtp.planned_los[bmtp.current_stay_num]
        yield self.env.timeout(los)

        # Go to next destination (which could be an exitflow)
        next_unit_name = bmtp.planned_route_stop[bmtp.current_stay_num + 1]
        self.out = bmtp_units[next_unit_name]
        if bmtp.current_stay_num == bmtp.route_length:
            # For ExitFlow object, no process needed
            self.out.put(bmtp)
        else:
            # Process for putting patient into next bed
            self.env.process(self.out.put(bmtp))

    # def get(self, unit):
    #     """ Release """
    #     return unit.get()

    def basic_stats_msg(self):
        """ Compute entries, exits, avg los and create summary message.


        Returns
        -------
        str
            Message with basic stats
        """

        if self.num_exits > 0:
            alos = self.tot_occ_time / self.num_exits
        else:
            alos = 0

        msg = "{:6}:\t Entries={}, Exits={}, Occ={}, ALOS={:4.2f}".\
            format(self.name, self.num_entries, self.num_exits,
                   self.unit.count, alos)
        return msg


class ExitFlow(object):
    """ Patients routed here when ready to exit.

        Patient objects put into a Store. Can be accessed later for stats
        and logs. A little worried about how big the Store will get.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        debug : boolean
            if true then patient details printed on arrival
    """

    def __init__(self, env, name, store_bmtp=True, debug=False):
        self.store = simpy.Store(env)
        self.env = env
        self.name = name
        self.store_bmtp = store_bmtp
        self.debug = debug
        self.num_exits = 0
        self.last_exit = 0.0

    def put(self, bmtp):

        if bmtp.bed_requests[bmtp.current_stay_num] is not None:
            previous_request = bmtp.bed_requests[bmtp.current_stay_num]
            previous_unit_name = bmtp.planned_route_stop[bmtp.current_stay_num]
            previous_unit = bmtp_units[previous_unit_name]
            previous_unit.unit.release(previous_request)
            previous_unit.num_exits += 1
            previous_unit.tot_occ_time += self.env.now - bmtp.entry_ts[
                bmtp.current_stay_num]
            bmtp.exit_ts[bmtp.current_stay_num - 1] = self.env.now

        self.last_exit = self.env.now
        self.num_exits += 1

        if self.debug:
            print(bmtp)

        # Store patient
        if self.store_bmtp:
            self.store.put(bmtp)

    def basic_stats_msg(self):
        """ Create summary message with basic stats on exits.


        Returns
        -------
        str
            Message with basic stats
        """

        msg = "{:6}:\t Exits={}, Last Exit={:10.2f}".format(self.name,
                                                            self.num_exits,
                                                            self.last_exit)

        return msg


class BMT_Patient(object):
    """ Models an BMT patient

        Parameters
        ----------
        arr_time : float
            Patient arrival time
        patient_id : int
            Unique patient id
        patient_type : int
            Patient type id (default 1). Currently just one patient type.
            In our prior research work we used a scheme with 11 patient types.
        arr_stream : int
            Arrival stream id (default 1). Currently there is just one arrival
            stream corresponding to the one patient generator class. In future,
            likely to be be multiple generators for generating random and
            scheduled arrivals.

    """

    def __init__(self, arr_time, patient_id, patient_type=1, arr_stream=1,
                 prng=RandomState(0)):
        self.arr_time = arr_time
        self.patient_id = patient_id
        self.patient_type = patient_type
        self.arr_stream = arr_stream

        self.name = 'Patient_{}'.format(patient_id)

        # Hard coding route, los and bed requests for now
        # Not sure how best to do routing related data structures.
        # Hack for now using combination of lists here, the out member
        # and the obunits dictionary.
        self.current_stay_num = 0
        self.route_length = 3

        self.planned_route_stop = []
        self.planned_route_stop.append(None)
        self.planned_route_stop.append("RFL")
        self.planned_route_stop.append("CN1")
        self.planned_route_stop.append("PTP")
        self.planned_route_stop.append("CN2")
        self.planned_route_stop.append("CTH")
        self.planned_route_stop.append("STM")
        self.planned_route_stop.append("TPT")
        self.planned_route_stop.append("AFO")
        self.planned_route_stop.append("MTC")
        self.planned_route_stop.append("AFG")
        self.planned_route_stop.append("LTF")
        self.planned_route_stop.append("EXIT")

        self.planned_los = []
        self.planned_los.append(None)
        self.planned_los.append(prng.exponential(MEAN_LOS_RFL))
        self.planned_los.append(prng.exponential(MEAN_LOS_CN1))
        self.planned_los.append(prng.exponential(MEAN_LOS_PTP))
        self.planned_los.append(prng.exponential(MEAN_LOS_CN2))
        self.planned_los.append(prng.exponential(MEAN_LOS_CTH))
        self.planned_los.append(prng.exponential(MEAN_LOS_STM))
        self.planned_los.append(prng.exponential(MEAN_LOS_TPT))
        self.planned_los.append(prng.exponential(MEAN_LOS_AFO))
        self.planned_los.append(prng.exponential(MEAN_LOS_MTC))
        self.planned_los.append(prng.exponential(MEAN_LOS_AFG))
        self.planned_los.append(prng.exponential(MEAN_LOS_LTF))


        # Since we have fixed route for now, just initialize full list to
        # hold bed requests
        self.bed_requests = [None for _ in range(self.route_length + 1)]
        self.request_entry_ts = [None for _ in range(self.route_length + 1)]
        self.entry_ts = [None for _ in range(self.route_length + 1)]
        self.exit_ts = [None for _ in range(self.route_length + 1)]

    def __repr__(self):
        return "patientid: {}, arr_stream: {}, time: {}". \
            format(self.patient_id, self.arr_stream, self.arr_time)


class BMT_PatientGenerator(object):
    """ Generates patients.

        Set the "out" member variable to resource at which patient generated.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        arr_rate : float
            Poisson arrival rate (expected number of arrivals per unit time)
        patient_type : int
            Patient type id (default 1). Currently just one patient type.
            In our prior research work we used a scheme with 11 patient types.
        arr_stream : int
            Arrival stream id (default 1). Currently there is just one arrival
            stream corresponding to the one patient generator class. In future,
            likely to be be multiple generators for generating random and
            scheduled arrivals
        initial_delay : float
            Starts generation after an initial delay. (default 0.0)
        stoptime : float
            Stops generation at the stoptime. (default Infinity)
        max_arrivals : int
            Stops generation after max_arrivals. (default Infinity)
        debug: bool
            If True, status message printed after
            each patient created. (default False)

    """

    def __init__(self, env, arr_rate, patient_type=1, arr_stream=1,
                 initial_delay=0,
                 stoptime=simpy.core.Infinity,
                 max_arrivals=simpy.core.Infinity, debug=False):

        self.env = env
        self.arr_rate = arr_rate
        self.patient_type = patient_type
        self.arr_stream = arr_stream
        self.initial_delay = initial_delay
        self.stoptime = stoptime
        self.max_arrivals = max_arrivals
        self.debug = debug
        self.out = None
        self.num_patients_created = 0

        self.prng = RandomState(RNG_SEED)

        self.action = env.process(
            self.run())  # starts the run() method as a SimPy process

    def run(self):
        """The patient generator.
        """
        # Delay for initial_delay
        yield self.env.timeout(self.initial_delay)
        # Main generator loop that terminates when stoptime reached
        while self.env.now < self.stoptime and \
                        self.num_patients_created < self.max_arrivals:
            # Delay until time for next arrival
            # Compute next interarrival time
            iat = self.prng.exponential(1.0 / self.arr_rate)
            yield self.env.timeout(iat)
            self.num_patients_created += 1
            # Create new patient
            bmtp = BMT_Patient(self.env.now, self.num_patients_created,
                            self.patient_type, self.arr_stream,
                            prng=self.prng)

            if self.debug:
                print("Patient {} created at {:.2f}.".format(
                    self.num_patients_created, self.env.now))

            # Set out member to BMT_unit object representing next destination
            self.out = bmtp_units[bmtp.planned_route_stop[1]]
            # Initiate process of requesting first bed in route
            self.env.process(self.out.put(bmtp))


# Initialize a simulation environment
simenv = simpy.Environment()

# Compute and display traffic intensities
rho_rfl = ARR_RATE * MEAN_LOS_RFL / CAPACITY_RFL
rho_cn1 = ARR_RATE * MEAN_LOS_CN1 / CAPACITY_CN1
rho_ptp = ARR_RATE * MEAN_LOS_PTP / CAPACITY_PTP
rho_cn2 = ARR_RATE * MEAN_LOS_CN2 / CAPACITY_CN2
rho_cth = ARR_RATE * MEAN_LOS_CTH / CAPACITY_CTH
rho_stm = ARR_RATE * MEAN_LOS_STM / CAPACITY_STM
rho_tpt = ARR_RATE * MEAN_LOS_TPT / CAPACITY_TPT
rho_afo = ARR_RATE * MEAN_LOS_AFO / CAPACITY_AFO
rho_mtc = ARR_RATE * MEAN_LOS_MTC / CAPACITY_MTC
rho_afg = ARR_RATE * MEAN_LOS_AFG / CAPACITY_AFG
rho_ltf = ARR_RATE * MEAN_LOS_LTF / CAPACITY_LTF

print("rho_obs: {:6.3f}\nrho_ldr: {:6.3f}\nrho_pp: {:6.3f}".format(rho_rfl,
                                                                   rho_cn1,
                                                                   rho_ptp,
                                                                   rho_cn2,
                                                                   rho_cth,
                                                                   rho_stm,
                                                                   rho_tpt,
                                                                   rho_afo,
                                                                   rho_mtc,
                                                                   rho_afg,
                                                                   rho_ltf))

# Create nursing units
rfl_unit = BMT_unit(simenv, 'RFL', CAPACITY_RFL, debug=False)
cn1_unit = BMT_unit(simenv, 'CN1', CAPACITY_CN1, debug=True)
ptp_unit = BMT_unit(simenv, 'PTP', CAPACITY_PTP, debug=False)
cn2_unit = BMT_unit(simenv, 'CN2', CAPACITY_CN2, debug=True)
cth_unit = BMT_unit(simenv, 'CTH', CAPACITY_CTH, debug=False)
stm_unit = BMT_unit(simenv, 'STM', CAPACITY_STM, debug=True)
tpt_unit = BMT_unit(simenv, 'TPT', CAPACITY_TPT, debug=False)
afo_unit = BMT_unit(simenv, 'AFO', CAPACITY_AFO, debug=True)
mtc_unit = BMT_unit(simenv, 'MTC', CAPACITY_MTC, debug=False)
ltf_unit = BMT_unit(simenv, 'LTF', CAPACITY_LTF, debug=True)



# Define system exit
exitflow = ExitFlow(simenv, 'EXIT', store_bmtp=False)

# Create dictionary of units keyed by name. This object can be passed along
# to other objects so that the units are accessible as patients "flow".
bmtp_units = {'RFL': rfl_unit, 'CN1': cn1_unit, 'PTP': ptp_unit,
'CN2': cn2_unit, 'CTH': cth_unit, 'STM': stm_unit,
'TPT': tpt_unit, 'AFO': afo_unit, 'MTC': mtc_unit,
'AFO': afo_unit, 'LTF': ltf_unit, 'EXIT': exitflow}

# Create a patient generator
bmtp_gen = BMT_PatientGenerator(simenv, ARR_RATE, debug=False)

# Routing logic
# Currently routing logic is hacked into the OBPatientGenerator
# and BMT
# Patient objects

# Run the simulation for a while
runtime = 1000
simenv.run(until=runtime)

# Patient generator stats
print("\nNum patients generated: {}\n".format(bmtp_gen.num_patients_created))

# Unit stats
print(rfl_unit.basic_stats_msg())
print(cn1_unit.basic_stats_msg())
print(ptp_unit.basic_stats_msg())
print(cn2_unit.basic_stats_msg())
print(cth_unit.basic_stats_msg())
print(stm_unit.basic_stats_msg())
print(tpt_unit.basic_stats_msg())
print(afo_unit.basic_stats_msg())
print(mtc_unit.basic_stats_msg())
print(ltf_unit.basic_stats_msg())

# System exit stats
print("\nNum patients exiting system: {}\n".format(exitflow.num_exits))
print("Last exit at: {:.2f}\n".format(exitflow.last_exit))
