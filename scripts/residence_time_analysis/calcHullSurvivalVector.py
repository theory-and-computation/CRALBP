import matplotlib.pyplot as plt
import sqlite3
from sksurv.nonparametric import kaplan_meier_estimator
import scienceplots
import matplotlib



def getContacts(start_resid, end_resid):
    connection = sqlite3.connect("/Users/danielsantos/Downloads/CRALBP/scripts/system_analysis_prep/RSCM_contacts.db")  # connect to your DB
    cursor = connection.cursor()  # get a cursor

    POPS_contacts = {contact: 0 for (contact, res) in cursor.execute(
                    f"""SELECT distinct frame, resid
                    FROM protein
                    WHERE resid BETWEEN {start_resid} and {end_resid}
                    AND lipid_restypes LIKE '%POPS%'
                    """)}

    POPC_contacts = {contact: 0 for (contact, res) in cursor.execute(
                    f"""SELECT distinct frame, resid
                    FROM protein
                    WHERE resid BETWEEN {start_resid} and {end_resid}
                    AND lipid_restypes LIKE '%POPC%'
                    """)}

    return POPC_contacts, POPS_contacts


def calcResidenceTimes(contact_frames, num_frames):
    residence_times = []
    stay = 1
    for frame in range(1, num_frames):
        if (frame - 1) in contact_frames:
            stay += 1
        elif (frame + 1) in contact_frames:
            residence_times.append(stay)
        else:
            stay = 1
    return residence_times


def calculateSurvivalVector(residence_times):
    dummy_event_titles = [True for x in residence_times]
    time_interval, survival_prob, conf_interval = kaplan_meier_estimator(event=dummy_event_titles, time_exit=residence_times, conf_type="log-log")
    return time_interval, survival_prob, conf_interval


def plotSurvival(popc_survival_vector, pops_survival_vector):

    fig, ax = plt.subplots(figsize = (15, 8))

    plt.plot(popc_survival_vector[0], popc_survival_vector[1], label='POPC', linewidth=2, color='#BDDBAA')
    plt.plot(pops_survival_vector[0], pops_survival_vector[1], label='POPS', linewidth=2, color='#C878AD')
    plt.fill_between(x=popc_survival_vector[0], y1=popc_survival_vector[2][0], y2=popc_survival_vector[2][1], alpha=0.25, color='#BDDBAA')
    plt.fill_between(x=pops_survival_vector[0], y1=pops_survival_vector[2][0], y2=pops_survival_vector[2][1], alpha=0.25, color='#C878AD')

    plt.gca().set_xlim(0, 1000)
    plt.gca().set_ylim(0, 1.0)

    ax.tick_params(which = 'major', length = 7)
    ax.tick_params(which = 'minor', direction = 'inout', length = 5, width = 1)
    ax.tick_params(which = 'both', bottom = True, top = True, left = True, right = True)

    figure = ax.get_figure()

    figure.savefig('hull_survival.pdf')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'Arial'
    plt.style.use(['science', 'nature', 'no-latex'])

    POPC_con_hull, POPS_con_hull = getContacts(252, 268)

    POPC_rt_hull = calcResidenceTimes(POPC_con_hull, 90000)
    POPS_rt_hull = calcResidenceTimes(POPS_con_hull, 90000)

    POPC_sv_hull = calculateSurvivalVector(POPC_rt_hull)
    POPS_sv_hull = calculateSurvivalVector(POPS_rt_hull)
    plotSurvival(POPC_sv_hull, POPS_sv_hull)
