import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

coeffs = [0.0020541, 0, -0.0207295, 0, 0.249792, 0.5, 0]

poly = np.poly1d(coeffs)

my_range = 5

# 生成样本数据
def relu(x):
    return np.maximum(0, x)

# x 的范围
x = np.linspace(-my_range, my_range, 400)
y = relu(x)
y_ = poly(x)

# 多项式拟合
degree = 3  # 多项式的阶数
p = Polynomial.fit(x, y, degree)
p = p.convert()  # 转换多项式以获取系数

# 生成拟合曲线的数据
x_fit = np.linspace(-my_range, my_range, 400)
y_fit = p(x_fit)

# 提取多项式系数
coeffs = p.coef
poly_str = " + ".join([f"{coeff:.4f}x^{i}" for i, coeff in enumerate(coeffs)])

# 打印多项式
print("The polynomial fit is:")
print(poly_str)

# 可视化
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='ReLU Function', linewidth=2)
plt.plot(x, y_, label='assumption', linewidth=2)
plt.plot(x_fit, y_fit, label=f'Polynomial Fit (degree={degree})', linestyle='--', linewidth=2)
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Polynomial Fit to ReLU Function')
plt.grid(True)
plt.show()


