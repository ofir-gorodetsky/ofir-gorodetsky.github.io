import numpy
import scipy.io
import matplotlib.pyplot as plt

def generate_all_graphs():
    graph_AP_variance()
    graphs_G()
    graphs_F()
    graph_AP_large_p()
    graph_short_interval_variance()

def z_measure(z, pp):
    #Computes the probability of a partition pp under the z-measure, but does not divide by (n! * rising factorial).
    dim_squared = pp.dimension()^2
    content_prod = prod([z + c[1] - c[0] for c in pp.cells()])^2
    return dim_squared * content_prod

def return_lambda1_dist(z = 0.5, n = 30):
    #Returns a vector of tuples (i,p), where:
    #i is of the form a/n, and p is the probability that lambda_1 = a, under the z-measure on Y_n.
    l_probs = [0 for i in xrange(n+1)]
    for par in Partitions(n):
        l_probs[par[0]] += z_measure(z,par)
    n_factorial = factorial(n)
    pochhammer = rising_factorial(z^2, n)
    v = [(i * 1.0 / n, l_probs[i] / (n_factorial * pochhammer) * 1.0) for i in xrange(len(l_probs))]
    return v

def graphs_F(n = 60):
    points_list = []
    points_deriv = []
    for z in [0.25, 0.5, 0.75, 1.5, 2.0, 2.5, 3.0, 3.5]:
        #The probability distribution of lambda_1/n under the n-th z-measure approximates the
        #probability distribution of alpha_1 under the z-measure.
        lambda1_dist = return_lambda1_dist(z, n=n)
        #Computing commulative distribution function of lambda_1 from the probability distribution function:
        lambda1_dist_cumsum = numpy.cumsum([t[1] for t in lambda1_dist]).tolist()
        lambda1_dist_CDF = [(lambda1_dist[i][0], lambda1_dist_cumsum[i]) for i in xrange(len(lambda1_dist))]

        points_list.append((z, lambda1_dist_CDF))
        #Computing discrete derivative of commulative distribution function:
        points_deriv.append((z,[(lambda1_dist_CDF[i][0], lambda1_dist_CDF[i+1][1]-lambda1_dist_CDF[i][1]) for i in xrange(len(lambda1_dist)-1)]))

    markers_smallz = ['_', '|', '+']
    markers_largez  =["x", "+", "D", "o", "d"]

    colors_smallz = [(0.0,0.0,0.5), (0.0,0.5,0.0), (0.5,0.0,0.0)]
    colors_largez = [(0.8,0.2,0.0), (0.8,0.5,0.0), (0.0,0.6,0.0), (0.0,0.5,0.8), (0.0,0.0,0.5)]

    #Graph of (approximation of) F_z for z = 0.25, z=0.5, z=0.75:
    fig = plt.figure()
    fig.patch.set_alpha(0.0)

    plt.scatter([s[0] for s in points_list[0][1]], [s[1] for s in points_list[0][1]], edgecolors='none', marker=markers_smallz[0], label='z=0.25', s=15, c=colors_smallz[0])
    plt.scatter([s[0] for s in points_list[1][1]], [s[1] for s in points_list[1][1]], edgecolors='none', marker=markers_smallz[1], label='z=0.5', s=15, c=colors_smallz[1])
    plt.scatter([s[0] for s in points_list[2][1]], [s[1] for s in points_list[2][1]], edgecolors='none', marker=markers_smallz[2], label='z=0.75', s=15, c=colors_smallz[2])
    plt.legend(loc='upper left', scatterpoints=1, fontsize=12, markerscale=1.5)

    axes = plt.gca()
    axes.set_xlim([0,1.01])
    axes.set_ylim([0,1.01])

    plt.savefig('Pz_n60_less1_shapes.png', dpi=100, bbox_inches='tight')
    plt.show()

    #Graph of (approximation of) F_z for z = 1.5, 2.0, 2.5, 3.0, 3.5:
    fig = plt.figure()
    fig.patch.set_alpha(0.0)

    plt.scatter([s[0] for s in points_list[7][1]], [s[1] for s in points_list[7][1]], edgecolors='none', marker=markers_largez[4], label='z=3.5', s=15, c=colors_largez[4])
    plt.scatter([s[0] for s in points_list[6][1]], [s[1] for s in points_list[6][1]], edgecolors='none', marker=markers_largez[3], label='z=3.0', s=15, c=colors_largez[3])
    plt.scatter([s[0] for s in points_list[5][1]], [s[1] for s in points_list[5][1]], edgecolors='none', marker=markers_largez[2], label='z=2.5', s=15, c=colors_largez[2])
    plt.scatter([s[0] for s in points_list[4][1]], [s[1] for s in points_list[4][1]], edgecolors='none', marker=markers_largez[1], label='z=2.0', s=15, c=colors_largez[1])
    plt.scatter([s[0] for s in points_list[3][1]], [s[1] for s in points_list[3][1]], edgecolors='none', marker=markers_largez[0], label='z=1.5', s=15, c=colors_largez[0])
    plt.legend(loc='upper left', scatterpoints=1, fontsize=12, markerscale=1.5)

    axes = plt.gca()
    axes.set_xlim([0,1.01])
    axes.set_ylim([0,1.01])

    plt.savefig('Pz_n60_bigger1_shapes.png', dpi=100, bbox_inches='tight')
    plt.show()

    #Graph of discrete derivative of (approximation of) F_z for z = 1.5, 2.0, 2.5, 3.0, 3.5:
    fig = plt.figure()
    fig.patch.set_alpha(0.0)

    plt.scatter([s[0] for s in points_deriv[7][1]], [s[1] for s in points_deriv[7][1]], edgecolors='none', marker=markers_largez[4], label='z=3.5', s=15, c=colors_largez[4])
    plt.scatter([s[0] for s in points_deriv[6][1]], [s[1] for s in points_deriv[6][1]], edgecolors='none', marker=markers_largez[3], label='z=3.0', s=15, c=colors_largez[3])
    plt.scatter([s[0] for s in points_deriv[5][1]], [s[1] for s in points_deriv[5][1]], edgecolors='none', marker=markers_largez[2], label='z=2.5', s=15, c=colors_largez[2])
    plt.scatter([s[0] for s in points_deriv[4][1]], [s[1] for s in points_deriv[4][1]], edgecolors='none', marker=markers_largez[1], label='z=2.0', s=15, c=colors_largez[1])
    plt.scatter([s[0] for s in points_deriv[3][1]], [s[1] for s in points_deriv[3][1]], edgecolors='none', marker=markers_largez[0], label='z=1.5', s=15, c=colors_largez[0])
    plt.legend(loc='upper left', scatterpoints=1, fontsize=12, markerscale=1.5)

    axes = plt.gca()
    axes.set_xlim([0,1.01])
    axes.set_ylim([0,0.13])

    plt.savefig('Pz_n60_bigger1_deriv_shapes.png', dpi=100, bbox_inches='tight')
    plt.show()

def return_G_dist(z = 0.5, n = 35):
    #Returns a vector (T(n;a))_{a=0,...n} a, where T(n;a) is defined in Theorem 1.1.
    lls = [[1 for i in xrange(n+2)]]
    for i in xrange(1, n+1):
        ll = map(lambda x: x[1], return_lambda1_dist(z,i))
        ll = numpy.cumsum(ll).tolist()
        lls.append(ll + [ll[-1] for i in xrange(n+2-len(ll))])
    #the list lls satisfies lls[i][a] = P_{Y_i}(\lambda_1 <= a), also for a>i.
    probs = [0] + [0 for i in xrange(n)]
    for i in xrange(n+1):
        j= n-i
        bin1 = falling_factorial(i+z^2-1,i)/factorial(i)
        bin2 = falling_factorial(j+z^2-1,j)/factorial(j)
        for a in xrange(1,len(probs)):
            probs[a] += bin1*bin2*lls[i][a-1]*lls[j][a]
    return probs

def graphs_G(n = 50):
    probs = return_G_dist(n=n)
    #We want to normalize T(n;a) by T(n;n). Asymptotically, T(n,n) ~ 1/sqrt(pi *n),
    #but in fact a simple closed form formula for T(n,n) exists:
    normalization_fact = binomial(n - 0.5, n) - binomial(n - 0.5, n)^2
    #The ratio T(n;a)/T(n;n) is important since it approximates G(a/n):
    G_approx = [(i * 1.0 / n, probs[i] * 1.0 / normalization_fact) for i in xrange(len(probs))]
    #Graph of (approximation of) G(s):
    g = Graphics()
    g += plot(points(G_approx, rgbcolor = (0.0,0.0,0.5), pointsize = 15))
    g.save('G_n_50.png')
    show(g)
    #Graph of (approximation of) derivative of G(s). Defined as the discrete derivatve of T(n;a)/T(n;n):
    g = Graphics()
    g += plot(points([(G_approx[i][0], G_approx[i+1][1]-G_approx[i][1]) for i in xrange(len(G_approx)-1)], rgbcolor = (0.0,0.0,0.5), pointsize = 15))
    g.save('G_n_50_deriv.png')
    show(g)

def graph_AP_variance(path = 'results.mat', n = 50):
    probs = return_G_dist(n=n)
    normalization_fact = binomial(n - 0.5, n) - binomial(n - 0.5, n)^2
    G_approx = [(i * 1.0 / n, probs[i] * 1.0 / normalization_fact) for i in xrange(len(probs))]

    mat = scipy.io.loadmat(path)
    x_axis = mat['x_axis']
    qcount = mat['qcount']
    K = 0.76422365 #The Landau-Ramanujan constant
    G_timesK = [(1-G_approx[i][0], G_approx[i][1] * K) for i in xrange(len(G_approx))]
    Empirical = [(x_axis[i], qcount[i]) for i in xrange(len(x_axis))]

    g = Graphics()
    g += plot(spline(G_timesK), 0, 1, rgbcolor = (0.0,0.0,0.5), legend_label = 'Prediction', legend_color = 'black')
    g += plot(points(Empirical, rgbcolor = (0.0,0.5,0.0), legend_label = 'Data', legend_color = 'black', pointsize = 3))
    show(g)
    g.save('Predicton_vs_Data_AP.png')

def graph_AP_large_p(path = 'results_large_p.mat'):
    mat = scipy.io.loadmat(path)
    x_axis = mat['x_axis']
    qcount = mat['qcount_large']

    K = 0.76422365 #The Landau-Ramanujan constant
    X = 9*10^8
    Prediction = [(delta, K + (-K^2 * X^delta + 1-X^(-delta))/(sqrt(math.log(X)))) for delta in [(math.log(2)/math.log(X))*j/100 for j in xrange(1,101)]]
    Empirical = [(x_axis[i], qcount[i]) for i in xrange(len(x_axis))]
    g = Graphics()
    g += plot(spline(Prediction), 0, math.log(2)/math.log(X), rgbcolor = (0.0,0.0,0.5), legend_label = 'Prediction', legend_color = 'black')
    g += plot(points(Empirical, rgbcolor = (0.0,0.5,0.0), legend_label = 'Data', legend_color = 'black', pointsize = 14))
    show(g)
    g.save('Predicton_vs_Data_AP_large_p.png')


def graph_short_interval_variance(path = 'results_short_interval.mat', n = 50):
    probs = return_G_dist(n=n)
    normalization_fact = binomial(n - 0.5, n) - binomial(n - 0.5, n)^2
    G_approx = [(i * 1.0 / n, probs[i] * 1.0 / normalization_fact) for i in xrange(len(probs))]

    mat = scipy.io.loadmat(path)
    x_axis = mat['x_axis']
    hcount = mat['hcount']
    K = 0.76422365 #The Landau-Ramanujan constant
    G_timesK = [(1-G_approx[i][0], G_approx[i][1] * K) for i in xrange(len(G_approx))]
    Empirical = [(x_axis[i], hcount[i]) for i in xrange(len(x_axis))]

    g = Graphics()
    g += plot(spline(G_timesK), 0, 1, rgbcolor = (0.0,0.0,0.5), legend_label = 'Prediction', legend_color = 'black')
    g += plot(points(Empirical, rgbcolor = (0.0,0.5,0.0), legend_label = 'Data', legend_color = 'black', pointsize = 3))
    show(g)
    g.save('Predicton_vs_Data_short_interval.png')