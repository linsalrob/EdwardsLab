import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats.distributions import norm
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity

data = [1.5] * 7 + [2.5] * 2 + [3.5] * 8 + [4.5] * 3 + [5.5] * 1 + [6.5] * 8
density = gaussian_kde(data)
xs = np.linspace(0, 8, 200)

x_grid = np.linspace(-4.5, 3.5, 1000)
x = np.concatenate([norm(-1, 1.).rvs(400), norm(1, 0.3).rvs(100)])

# print(x)

# print("X Grid:" + str(len(x_grid)))
# print("X: " + str(len(x)))

if (0):
    density.covariance_factor = lambda: .25
    density._compute_covariance()
    plt.plot(xs, density(xs))
    # plt.plot(data)
    plt.show()


bandwidth=0.2
kde_skl = KernelDensity(bandwidth=bandwidth)
kde_skl.fit(x[:, np.newaxis])
# score_samples() returns the log-likelihood of the samples
log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
print(log_pdf[0:10])
density = np.exp(log_pdf)


plt.plot(x_grid, density, color='blue', alpha=0.5, lw=3)
#plt.show()


grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.linspace(0.1, 1.0, 30)},
                    cv=20)  # 20-fold cross-validation
# print(x[:, None])
grid.fit(x[:, None])
# print grid.best_params_