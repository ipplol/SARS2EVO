# author
# https://gist.github.com/ruoyu0088/70effade57483355bbd18b31dc370f2a
import pylab as pl
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt


def make_test_data(seg_count, point_count):
    x = np.random.uniform(2, 10, seg_count)
    x = np.cumsum(x)
    x *= 10 / x.max()
    y = np.cumsum(np.random.uniform(-1, 1, seg_count))
    X = np.random.uniform(0, 10, point_count)
    Y = np.interp(X, x, y) + np.random.normal(0, 0.05, point_count)
    return X, Y


def segments_fit(X, Y, count):
    xmin = X.min()
    xmax = X.max()

    seg = np.full(count - 1, (xmax - xmin) / count)
    # print('seg:', seg)

    px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
    py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.01].mean() for x in px_init])
    # print('px_init:', px_init, 'py_init', py_init)

    def func(p):
        seg = p[:count - 1]
        py = p[count - 1:]
        px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        return px, py

    def err(p):
        px, py = func(p)
        Y2 = np.interp(X, px, py)
        return np.mean((Y - Y2)**2)

    r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    # print('r:',r)
    # print('r.xcol:', r.xcol)
    return func(r.x)

# OPLRA (Optimal Piecewise Linear Regression Analysis)
# To achieve this result, information criteria such as the Akaike and
# the Bayesian are used that reward predictive accuracy and penalise model complexity.
# optimizing AIC and/ or BIC
# https://discovery.ucl.ac.uk/id/eprint/10070516/1/AIC_BIC_Paper.pdf
def auto_segments_fit(X, Y, maxcount):
    xmin = X.min()
    xmax = X.max()

    n = len(X)

    AIC_ = float('inf')
    BIC_ = float('inf')
    r_ = None

    for count in range(1, maxcount + 1):

        seg = np.full(count - 1, (xmax - xmin) / count)

        px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.1].mean() for x in px_init])

        def func(p):
            seg = p[:count - 1]
            py = p[count - 1:]
            px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
            return px, py

        def err(p):  # This is RSS / n
            px, py = func(p)
            Y2 = np.interp(X, px, py)
            return np.mean((Y - Y2) ** 2)

        r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')

        # Compute AIC/ BIC.
        AIC = n * np.log10(err(r.x)) + 4 * count
        BIC = n * np.log10(err(r.x)) + 2 * count * np.log(n)

        if ((BIC < BIC_) & (AIC < AIC_)):  # Continue adding complexity.
            r_ = r
            AIC_ = AIC
            BIC_ = BIC
        else:  # Stop.
            count = count - 1
            break

    return func(r_.x)  ## Return the last (n-1)

def calculate_slope(x, y):
    slopes = []
    for i in range(x.shape[0] - 1):
        adj_point_x= np.r_[x[i], x[i+1]]
        adj_point_y = np.r_[y[i], y[i+1]]
        slope, intercept = np.polyfit(adj_point_x, adj_point_y, 1)
        slopes.append(slope)
    # print(slopes)
    return slopes


def fit_curve():

    def func(x, k, b):
        return k*x + b

    xdata = np.linspace(0, 10, 80)
    y = func(xdata, 2, 10)
    rng = np.random.default_rng()
    y_noise = 0.5 * rng.normal(size=xdata.size)
    ydata = y + y_noise
    plt.plot(xdata, ydata, 'ob', label='data')

    popt, pcov = optimize.curve_fit(func, xdata, ydata)

    plt.plot(xdata, func(xdata, *popt), 'r-',
             label='fit: a=%5.3f, b=%5.3f '% tuple(popt))

    popt, pcov = optimize.curve_fit(func, xdata, ydata, bounds=(0, [np.inf, 0.1]))

    plt.plot(xdata, func(xdata, *popt), 'g--',
             label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

    plt.xlabel('xcol')
    plt.ylabel('ycol')
    plt.legend()
    plt.show()

def np_polyfit(x, y, plot=False, ax=None):
    coefficient = np.polyfit(x, y, 1)
    p = np.poly1d(coefficient)

    if ax:
        handle = ax
    else:
        handle = mpl.pyplot

    if plot:
        xp = np.linspace(x.min(), x.max(), 100)
        handle.plot(xp, p(xp), '--', c='grey', lw=3)

    return coefficient[0]

def np_polyfit_derivative(x, y, plot=False, order=1):
    coefficient = np.polyfit(x, y, order)
    p = np.poly1d(coefficient)
    derivative = p.deriv()

    if plot:
        fig, axs = plt.subplots(1,2, figsize=(10, 6))
        xp = np.linspace(x.min(), x.max(), 100)
        axs[0].plot(xp, p(xp), '--', c='r')
        axs[0].plot(x, y, 'ob', label='data')
        axs[0].plot(xp, derivative(xp), '--', c='brown')
        axs[1].plot(xp, derivative(xp), '--', c='brown')
        plt.show()

    return derivative



if __name__ == '__main__':

    X, Y = make_test_data(8, 2000)
    px, py = segments_fit(X, Y, 8)
    slopes = calculate_slope(px, py)

    pl.plot(X, Y, ".")
    pl.plot(px, py, "-or")
    texts = [pl.text((px[i]+px[i+1])/2, (py[i]+py[i+1])/2, f'{slopes[i]:.2f}') for i in range(len(slopes))]
    pl.show()

    # X = np.random.uniform(0, 10, 2000)
    # Y = np.sin(X) + np.random.normal(0, 0.05, X.shape)
    # px, py = segments_fit(X, Y, 8)
    # pl.plot(X, Y, ".")
    # pl.plot(px, py, "-or")
    # pl.show()

    # path_results = Path('results/cov/mutation')
    # my_nextstrain_metadata = pd.read_csv(path_results/'my_nextstrain_metadata.csv')
    #
    # X = my_nextstrain_metadata['days'].values
    # Y = my_nextstrain_metadata['mutation_count'].values
    #
    # px, py = segments_fit(X, Y, 8)
    # pl.plot(X, Y, ".")
    # pl.plot(px, py, "-or")
    # pl.show()
    #
    # px, py = auto_segments_fit(X, Y, 10)
    # pl.plot(X, Y, ".")
    # pl.plot(px, py, "-or")
    # pl.show()
    #
    # for i in range(1,11):
    #     px, py = segments_fit(X, Y, i)
    #     pl.plot(X, Y, ".")
    #     pl.plot(px, py, "-or", linewidth=3, markersize=8)
    #     pl.title(f'segment = {i}')
    #     pl.show()

    # fit_curve()

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    x = np.linspace(0, 10, 80)
    y = func(x, 2.5, 1.3, 0.5)
    deriv = np_polyfit_derivative(x, y, plot=True, order=6)










