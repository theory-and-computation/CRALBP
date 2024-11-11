import matplotlib.pyplot as plt
import sqlite3
import numpy as npy
import scienceplots
import matplotlib
from scipy.stats import sem



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


def calculateMeanResidenceTime(residence_times):
    return npy.mean(residence_times)


def plotMean(POPC_mean, POPS_mean, POPC_stderr, POPS_stderr):
    fig, ax = plt.subplots(figsize = (15, 8))

    ax.spines['bottom'].set_linewidth(2)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_linewidth(0)
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_linewidth(0)
    ax.spines['right'].set_color('black')

    plt.bar(x=['POPC'], height=[POPC_mean], yerr=POPC_stderr, color='#BDDBAA', capsize=10, linewidth=2)
    plt.bar(x=['POPS'], height=[POPS_mean], yerr=POPS_stderr, color='#C878AD', capsize=10, linewidth=2)
    plt.gca().set_ylim(0, 350)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)

    ax.tick_params(which = 'both', bottom = True, top = False, left = True, right = False)

    figure = ax.get_figure()
    figure.savefig('hull_mean.pdf')
    plt.show()


if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'Arial'
    plt.style.use(['science', 'nature', 'no-latex'])

    POPC_con_hull, POPS_con_hull = getContacts(252, 268)

    POPC_rt_hull = calcResidenceTimes(POPC_con_hull, 90000)
    POPS_rt_hull = calcResidenceTimes(POPS_con_hull, 90000)


    POPC_mrt = calculateMeanResidenceTime(POPC_rt_hull)
    POPS_mrt = calculateMeanResidenceTime(POPS_rt_hull)

    plotMean(POPC_mrt, POPS_mrt, sem(POPC_rt_hull), sem(POPS_rt_hull))
