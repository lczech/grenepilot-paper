import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(42)


# Population level quantities
def population_pi_within(freqs1, freqs2):
    return 1 - 0.5 * ((freqs1**2).sum() + (freqs2**2).sum())


def population_pi_between(freqs1, freqs2):
    return 1 - freqs1.dot(freqs2)


def population_pi_total(freqs1, freqs2):
    return 0.5 * (population_pi_within(freqs1, freqs2)
                  + population_pi_between(freqs1, freqs2))


def population_fst_nei(freqs1, freqs2):
    return 1 - (population_pi_within(freqs1, freqs2)
                / population_pi_total(freqs1, freqs2))


def population_fst_hudson(freqs1, freqs2):
    return 1 - (population_pi_within(freqs1, freqs2)
                / population_pi_between(freqs1, freqs2))


# Simulation functions
def simulate(freqs1, freqs2, n1, n2, c1, c2):
    pool_freq_1 = np.random.multinomial(n1, freqs1) / n1
    pool_freq_2 = np.random.multinomial(n2, freqs2) / n2
    read_freq_1 = np.random.multinomial(c1, pool_freq_1) / c1
    read_freq_2 = np.random.multinomial(c2, pool_freq_2) / c2
    return read_freq_1, read_freq_2


def batch_simulate_random(reps, n1, n2, c1, c2):
    truth_1 = []
    truth_2 = []
    sample_1 = []
    sample_2 = []
    for _ in range(reps):
        f1 = np.random.dirichlet([4, 1, 1, 1])
        f2 = np.random.dirichlet([4, 1, 1, 1])
        s1, s2 = simulate(f1, f2, n1, n2, c1, c2)
        truth_1.append(f1)
        truth_2.append(f2)
        sample_1.append(s1)
        sample_2.append(s2)
    truth_1 = np.array(truth_1)
    truth_2 = np.array(truth_2)
    sample_1 = np.array(sample_1)
    sample_2 = np.array(sample_2)
    return truth_1, truth_2, sample_1, sample_2


def batch_simulate_fixed(reps, freqs1, freqs2, n1, n2, c1, c2):
    sample_1 = []
    sample_2 = []
    for _ in range(reps):
        s1, s2 = simulate(freqs1, freqs2, n1, n2, c1, c2)
        sample_1.append(s1)
        sample_2.append(s2)
    return np.array(sample_1), np.array(sample_2)


# Unbiased estimators of pi
def unbiased_pi_within(f1, f2, n1, n2, c1, c2):
    pi_1 = n1/(n1-1)*c1/(c1-1)*(1 - (f1**2).sum())
    pi_2 = n2/(n2-1)*c2/(c2-1)*(1 - (f2**2).sum())
    return 0.5*(pi_1 + pi_2)


def unbiased_pi_between(f1, f2, n1, n2, c1, c2):
    return 1 - f1.dot(f2)


def unbiased_pi_total(f1, f2, n1, n2, c1, c2):
    return 0.5*(unbiased_pi_within(f1, f2, n1, n2, c1, c2)
                + unbiased_pi_between(f1, f2, n1, n2, c1, c2))


# Karlsson estimators of pi
def karlsson_pi_within(f1, f2, n1, n2, c1, c2):
    pi_1 = c1/(c1-1)*(1 - (f1**2).sum())
    pi_2 = c2/(c2-1)*(1 - (f2**2).sum())
    pi_within = 0.5*(pi_1 + pi_2)
    return pi_within


def karlsson_pi_between(f1, f2, n1, n2, c1, c2):
    pi_between = 1 - f1.dot(f2)
    return pi_between


# "Conventional" estimators of pi
def conventional_pi_within(f1, f2, n1, n2, c1, c2):
    pi_1 = n1/(n1-1)*c1/(c1-1)*(1 - (f1**2).sum())
    pi_2 = n2/(n2-1)*c2/(c2-1)*(1 - (f2**2).sum())
    return 0.5*(pi_1+pi_2)


def conventional_pi_total(f1, f2, n1, n2, c1, c2):
    c_min = np.min([c1, c2])
    n_min = np.min([n1, n2])
    f_t = 0.5 * (f1 + f2)
    pi_t = n_min/(n_min-1)*c_min/(c_min-1)*(1 - (f_t**2).sum())
    return pi_t


# All of the Fsts
def karlsson_hudson(f1, f2, n1, n2, c1, c2):
    pi_within = karlsson_pi_within(f1, f2, n1, n2, c1, c2)
    pi_between = karlsson_pi_between(f1, f2, n1, n2, c1, c2)
    return 1 - pi_within / pi_between


def our_hudson(f1, f2, n1, n2, c1, c2):
    return 1 - (unbiased_pi_within(f1, f2, n1, n2, c1, c2)
                / unbiased_pi_between(f1, f2, n1, n2, c1, c2))


def conventional_nei(f1, f2, n1, n2, c1, c2):
    pi_w = conventional_pi_within(f1, f2, n1, n2, c1, c2)
    pi_t = conventional_pi_total(f1, f2, n1, n2, c1, c2)
    return 1 - pi_w / pi_t


def our_nei(f1, f2, n1, n2, c1, c2):
    return 1 - (unbiased_pi_within(f1, f2, n1, n2, c1, c2)
                / unbiased_pi_total(f1, f2, n1, n2, c1, c2))


def analyze_and_plot(n, c):
    num_reps = 100001
    window_size = 100
    print('Simulating n={}, c={}'.format(n, c))
    # simulate random data
    ts1, ts2, ss1, ss2 = batch_simulate_random(num_reps,
                                               n1=n, n2=n, c1=c, c2=c)
    hudson_true = []
    nei_true = []
    karlsson_ests = []
    our_hudson_ests = []
    conventional_ests = []
    our_nei_ests = []
    true_pi_w = []
    true_pi_b = []
    true_pi_t = []
    our_pi_w = []
    our_pi_b = []
    our_pi_t = []
    conventional_pi_w = []
    conventional_pi_t = []
    karlsson_pi_w = []
    karlsson_pi_b = []

    curr_true_pi_w = []
    curr_true_pi_b = []
    curr_true_pi_t = []
    curr_our_pi_w = []
    curr_our_pi_b = []
    curr_our_pi_t = []
    curr_conventional_pi_w = []
    curr_conventional_pi_t = []
    curr_karlsson_pi_w = []
    curr_karlsson_pi_b = []
    for idx, (t1, t2, s1, s2) in enumerate(zip(ts1, ts2, ss1, ss2)):
        if idx > 0 and (idx % window_size) == 0:
            true_pi_w.append(np.mean(curr_true_pi_w))
            curr_true_pi_w = []

            true_pi_b.append(np.mean(curr_true_pi_b))
            curr_true_pi_b = []

            true_pi_t.append(np.mean(curr_true_pi_t))
            curr_true_pi_t = []

            our_pi_w.append(np.mean(curr_our_pi_w))
            curr_our_pi_w = []

            our_pi_b.append(np.mean(curr_our_pi_b))
            curr_our_pi_b = []

            our_pi_t.append(np.mean(curr_our_pi_t))
            curr_our_pi_t = []

            conventional_pi_w.append(np.mean(curr_conventional_pi_w))
            curr_conventional_pi_w = []

            conventional_pi_t.append(np.mean(curr_conventional_pi_t))
            curr_conventional_pi_t = []

            karlsson_pi_w.append(np.mean(curr_karlsson_pi_w))
            curr_karlsson_pi_w = []

            karlsson_pi_b.append(np.mean(curr_karlsson_pi_b))
            curr_karlsson_pi_b = []

        curr_true_pi_w.append(population_pi_within(t1, t2))
        curr_true_pi_b.append(population_pi_between(t1, t2))
        curr_true_pi_t.append(population_pi_total(t1, t2))

        curr_our_pi_w.append(unbiased_pi_within(s1, s2, n, n, c, c))
        curr_our_pi_b.append(unbiased_pi_between(s1, s2, n, n, c, c))
        curr_our_pi_t.append(unbiased_pi_total(s1, s2, n, n, c, c))

        curr_conventional_pi_w.append(
            conventional_pi_within(s1, s2, n, n, c, c)
        )
        curr_conventional_pi_t.append(
            conventional_pi_total(s1, s2, n, n, c, c)
        )

        curr_karlsson_pi_w.append(karlsson_pi_within(s1, s2, n, n, c, c))
        curr_karlsson_pi_b.append(karlsson_pi_between(s1, s2, n, n, c, c))

        hudson_true.append(population_fst_hudson(t1, t2))
        nei_true.append(population_fst_nei(t1, t2))
        karlsson_ests.append(karlsson_hudson(s1, s2, n, n, c, c))
        our_hudson_ests.append(our_hudson(s1, s2, n, n, c, c))
        conventional_ests.append(conventional_nei(s1, s2, n, n, c, c))
        our_nei_ests.append(our_nei(s1, s2, n, n, c, c))
    results = pd.DataFrame()
    results['True Hudson Fst'] = hudson_true
    results['True Nei Fst'] = nei_true
    results['Karlsson Estimator'] = karlsson_ests
    results['Our Hudson Fst Estimator'] = our_hudson_ests
    results['Conventional Estimator'] = conventional_ests
    results['Our Nei Fst Estimator'] = our_nei_ests

    results_windowed = pd.DataFrame()
    results_windowed['True Hudson Fst'] = [
        1 - pw/pb for pw, pb in zip(true_pi_w, true_pi_b)
    ]
    results_windowed['True Nei Fst'] = [
        1 - pw/pt for pw, pt in zip(true_pi_w, true_pi_t)
    ]
    results_windowed['Karlsson Estimator'] = [
        1 - pw/pb for pw, pb in zip(karlsson_pi_w, karlsson_pi_b)
    ]
    results_windowed['Our Hudson Fst Estimator'] = [
        1 - pw/pb for pw, pb in zip(our_pi_w, our_pi_b)
    ]
    results_windowed['Conventional Estimator'] = [
        1 - pw/pt for pw, pt in zip(conventional_pi_w, conventional_pi_t)
    ]
    results_windowed['Our Nei Fst Estimator'] = [
        1 - pw/pt for pw, pt in zip(our_pi_w, our_pi_t)
    ]

    for res, name in zip([results, results_windowed], ['SNP', 'window']):
        sns.regplot(data=res, x='True Hudson Fst', y='Karlsson Estimator',
                    scatter_kws={'s': 1},
                    line_kws={'linewidth': 4,
                              'ls': '--'})
        plt.plot([min(res['True Hudson Fst']), max(res['True Hudson Fst'])],
                 [min(res['True Hudson Fst']), max(res['True Hudson Fst'])],
                 c='red')
        plt.title('Karlsson Estimator at the {} level'.format(name))
        plt.savefig(
            'n_{}_c_{}_random_karlsson_vs_hudson_{}.png'.format(n, c, name)
        )
        plt.clf()

        sns.regplot(data=res, x='True Hudson Fst',
                    y='Our Hudson Fst Estimator',
                    scatter_kws={'s': 1},
                    line_kws={'linewidth': 4,
                              'ls': '--'})
        plt.plot([min(res['True Hudson Fst']), max(res['True Hudson Fst'])],
                 [min(res['True Hudson Fst']), max(res['True Hudson Fst'])],
                 c='red')
        plt.title('Our Hudson Fst Estimator at the {} level'.format(name))
        plt.savefig(
            'n_{}_c_{}_random_unbiased_vs_hudson_{}.png'.format(n, c, name)
        )
        plt.clf()

        sns.regplot(data=res, x='True Nei Fst', y='Conventional Estimator',
                    scatter_kws={'s': 1},
                    line_kws={'linewidth': 4,
                              'ls': '--'})
        plt.plot([min(res['True Nei Fst']), max(res['True Nei Fst'])],
                 [min(res['True Nei Fst']), max(res['True Nei Fst'])],
                 c='red')
        plt.title('Conventional Estimator at the {} level'.format(name))
        plt.savefig(
            'n_{}_c_{}_random_conventional_vs_nei_{}.png'.format(n, c, name)
        )
        plt.clf()
        sns.regplot(data=res, x='True Nei Fst', y='Our Nei Fst Estimator',
                    scatter_kws={'s': 1},
                    line_kws={'linewidth': 4,
                              'ls': '--'})
        plt.plot([min(res['True Nei Fst']), max(res['True Nei Fst'])],
                 [min(res['True Nei Fst']), max(res['True Nei Fst'])],
                 c='red')
        plt.title('Our Nei Fst Estimator at the {} level'.format(name))
        plt.savefig(
            'n_{}_c_{}_random_unbiased_vs_nei_{}.png'.format(n, c, name)
        )
        plt.clf()

    # Onto a fixed example
    f1 = np.array([0.2, 0.8, 0., 0.])
    f2 = np.array([0.3, 0.7, 0., 0.])
    hudson_true = population_fst_hudson(f1, f2)
    nei_true = population_fst_nei(f1, f2)
    karlsson_ests = []
    our_hudson_ests = []
    conventional_ests = []
    our_nei_ests = []

    true_pi_w = []
    true_pi_b = []
    true_pi_t = []
    our_pi_w = []
    our_pi_b = []
    our_pi_t = []
    conventional_pi_w = []
    conventional_pi_t = []
    karlsson_pi_w = []
    karlsson_pi_b = []

    curr_true_pi_w = []
    curr_true_pi_b = []
    curr_true_pi_t = []
    curr_our_pi_w = []
    curr_our_pi_b = []
    curr_our_pi_t = []
    curr_conventional_pi_w = []
    curr_conventional_pi_t = []
    curr_karlsson_pi_w = []
    curr_karlsson_pi_b = []

    ss1, ss2 = batch_simulate_fixed(num_reps, f1, f2, n, n, c, c)

    for idx, (s1, s2) in enumerate(zip(ss1, ss2)):
        if idx > 0 and idx % window_size == 0:
            true_pi_w.append(np.mean(curr_true_pi_w))
            curr_true_pi_w = []

            true_pi_b.append(np.mean(curr_true_pi_b))
            curr_true_pi_b = []

            true_pi_t.append(np.mean(curr_true_pi_t))
            curr_true_pi_t = []

            our_pi_w.append(np.mean(curr_our_pi_w))
            curr_our_pi_w = []

            our_pi_b.append(np.mean(curr_our_pi_b))
            curr_our_pi_b = []

            our_pi_t.append(np.mean(curr_our_pi_t))
            curr_our_pi_t = []

            conventional_pi_w.append(np.mean(curr_conventional_pi_w))
            curr_conventional_pi_w = []

            conventional_pi_t.append(np.mean(curr_conventional_pi_t))
            curr_conventional_pi_t = []

            karlsson_pi_w.append(np.mean(curr_karlsson_pi_w))
            curr_karlsson_pi_w = []

            karlsson_pi_b.append(np.mean(curr_karlsson_pi_b))
            curr_karlsson_pi_b = []

        karlsson_ests.append(karlsson_hudson(s1, s2, n, n, c, c))
        our_hudson_ests.append(our_hudson(s1, s2, n, n, c, c))
        conventional_ests.append(conventional_nei(s1, s2, n, n, c, c))
        our_nei_ests.append(our_nei(s1, s2, n, n, c, c))

        curr_our_pi_w.append(unbiased_pi_within(s1, s2, n, n, c, c))
        curr_our_pi_b.append(unbiased_pi_between(s1, s2, n, n, c, c))
        curr_our_pi_t.append(unbiased_pi_total(s1, s2, n, n, c, c))

        curr_conventional_pi_w.append(
            conventional_pi_within(s1, s2, n, n, c, c)
        )
        curr_conventional_pi_t.append(
            conventional_pi_total(s1, s2, n, n, c, c)
        )

        curr_karlsson_pi_w.append(karlsson_pi_within(s1, s2, n, n, c, c))
        curr_karlsson_pi_b.append(karlsson_pi_between(s1, s2, n, n, c, c))

    # SNP level
    plt.hist(karlsson_ests, bins=100)
    plt.plot([hudson_true]*2, [0, num_reps/20], c='red', linewidth=4)
    plt.plot([np.nanmean(karlsson_ests)]*2, [0, num_reps/20],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Karlsson Estimator')
    plt.title('Karlsson Estimator at the SNP level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_karlsson_vs_hudson_SNP.png'.format(n, c))
    plt.clf()

    plt.hist(our_hudson_ests, bins=100)
    plt.plot([hudson_true]*2, [0, num_reps/20], c='red', linewidth=4)
    plt.plot([np.nanmean(our_hudson_ests)]*2, [0, num_reps/20],
             c='black', ls='--', linewidth=4)

    plt.xlabel('Our Hudson Fst Estimator')
    plt.title('Our Hudson Fst Estimator at the SNP level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_unbiased_vs_hudson_SNP.png'.format(n, c))
    plt.clf()

    plt.hist(conventional_ests, bins=100)
    plt.plot([nei_true]*2, [0, num_reps/20], c='red', linewidth=4)
    plt.plot([np.nanmean(conventional_ests)]*2, [0, num_reps/20],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Conventional Estimator')
    plt.title('Conventional Estimator at the SNP level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_conventional_vs_nei_SNP.png'.format(n, c))
    plt.clf()

    plt.hist(our_nei_ests, bins=100)
    plt.plot([nei_true]*2, [0, num_reps/20], c='red', linewidth=4)
    plt.plot([np.nanmean(our_nei_ests)]*2, [0, num_reps/20],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Our Nei Fst Estimator')
    plt.title('Our Nei Fst Estimator at the SNP level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_unbiased_vs_nei_SNP.png'.format(n, c))
    plt.clf()

    # Window Level
    window_karlsson_ests = [
        1 - pw/pb for pw, pb in zip(karlsson_pi_w, karlsson_pi_b)
    ]
    plt.hist(window_karlsson_ests, bins=100)
    plt.plot([hudson_true]*2, [0, num_reps/20/window_size],
             c='red', linewidth=4)
    plt.plot([np.nanmean(window_karlsson_ests)]*2,
             [0, num_reps/20/window_size],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Karlsson Estimator')
    plt.title('Karlsson Estimator at the window level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_karlsson_vs_hudson_window.png'.format(n, c))
    plt.clf()

    window_our_hudson_ests = [
        1 - pw/pb for pw, pb in zip(our_pi_w, our_pi_b)
    ]
    plt.hist(window_our_hudson_ests, bins=100)
    plt.plot([hudson_true]*2, [0, num_reps/20/window_size],
             c='red', linewidth=4)
    plt.plot([np.nanmean(window_our_hudson_ests)]*2,
             [0, num_reps/20/window_size],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Our Hudson Fst Estimator')
    plt.title('Our Hudson Fst Estimator at the window level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_unbiased_vs_hudson_window.png'.format(n, c))
    plt.clf()

    window_conventional_ests = [
        1 - pw/pt for pw, pt in zip(conventional_pi_w, conventional_pi_t)
    ]
    plt.hist(window_conventional_ests, bins=100)
    plt.plot([nei_true]*2, [0, num_reps/20/window_size], c='red', linewidth=4)
    plt.plot([np.nanmean(window_conventional_ests)]*2,
             [0, num_reps/20/window_size],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Conventional Estimator')
    plt.title('Conventional Estimator at the window level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_conventional_vs_nei_window.png'.format(n, c))
    plt.clf()

    window_our_nei_ests = [
        1 - pw/pt for pw, pt in zip(our_pi_w, our_pi_t)
    ]
    plt.hist(window_our_nei_ests, bins=100)
    plt.plot([nei_true]*2, [0, num_reps/20/window_size], c='red', linewidth=4)
    plt.plot([np.nanmean(window_our_nei_ests)]*2,
             [0, num_reps/20/window_size],
             c='black', ls='--', linewidth=4)
    plt.xlabel('Our Nei Fst Estimator')
    plt.title('Our Nei Fst Estimator at the window level')
    plt.legend(['True Value', 'Mean of Estimates'])
    plt.savefig('n_{}_c_{}_fixed_unbiased_vs_nei_window.png'.format(n, c))
    plt.clf()


def _main():
    # large n, large c
    analyze_and_plot(n=1000, c=1000)

    # large n, small c
    analyze_and_plot(n=1000, c=10)

    # small n, large c
    analyze_and_plot(n=10, c=1000)

    # small n, small c
    analyze_and_plot(n=10, c=10)


if __name__ == '__main__':
    _main()
