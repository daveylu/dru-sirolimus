import matplotlib.pyplot as plt

x = [i for i in range(10)]
y = [i * 2 for i in range(10)]

plt.figure(figsize = (2560 / 96, 1440 / 96))
plt.scatter(x, y)
plt.show()