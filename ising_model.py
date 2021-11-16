#Встроенные
import argparse
import math
#Не встроенные
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

#====================Парсер аргументов из командной строки====================#
def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--gr-name', type=str, default='gr.png', help='Название графика (расширение обязательно)')
	parser.add_argument('-l', type=int, default=16, help='Размер квадратного поля: lxl (целое число)')

	args = parser.parse_args()
	cfg = vars(args)
	args = argparse.Namespace(**cfg)
	print(args, '\n')

	return args
#==============================Инициализатор поля==============================#
def init(L):
		state = 2 * np.random.randint(2, size=(L,L)) - 1
		return state

def magnetization(config):
	Mag = np.sum(config)
	return Mag
#=====================Алгоритм Метрополиса + метод Монте-Карло=================#
def MC_step(config, beta):
	L = len(config)
	for i in range(L):
		for j in range(L):
			a = np.random.randint(0, L)
			b = np.random.randint(0, L)
			sigma =  config[a, b]
			neighbors = config[(a+1)%L, b] + config[a, (b+1)%L] + config[(a-1)%L, b] + config[a, (b-1)%L]
			del_E = 2*sigma*neighbors
			if del_E < 0:
				sigma *= -1
			elif rand(1) < np.exp(-del_E*beta):
				sigma *= -1
			config[a, b] = sigma
	return config

def calcul_energy_mag_C_X(config, L, eqSteps=1, err_runs=1):

	nt      = 100  
	mcSteps = 1000
	T_c = 2/math.log(1 + math.sqrt(2))
	T = np.linspace(1., 7., nt); 
	M = np.zeros(nt)
	M_theoric = np.zeros(nt)
	n1 = 1.0/(mcSteps*L*L)
	n2 = 1.0/(mcSteps*mcSteps*L*L)
		
	Magnetizations = []
	for t in range(nt):
		beta = 1./T[t]
		for i in range(eqSteps):
			MC_step(config, beta)
		Mz = []

		for j in range(err_runs):
			M = M_squared = 0
			for i in range(mcSteps):
				MC_step(config, beta)           
				mag = abs(magnetization(config))
				M += mag
				M_squared += mag**2

			M_mean = M/mcSteps
			M_squared_mean = M_squared/mcSteps
			Magnetization = M_mean/L**2
			Mz.append(Magnetization)
		Magnetization = np.mean(Mz)
		Magnetizations.append(Magnetization)
		
		if T[t] - T_c < 0:
			M_theoric[t] = pow(1 - pow(np.sinh(2*beta), -4),1/8)

	return T, Magnetizations, M_theoric
#================Отрисовка гафика=================#	
def plot(args, T, Magnetizations, M_theoric):
	fig, ax = plt.subplots(nrows=1, ncols=1)
	ax.set_title("Зависимость намагниченности от температуры в модели Изинга", size=10)
	ax.set_xlabel('Температура', size=10)
	ax.set_ylabel('Намагниченность', size=10)
	ax.scatter(T, Magnetizations, color='g', marker='.', label='Размер поля: '+str(L)+'x'+str(L))
	ax.plot(T, M_theoric, color='red', label='Теория')
	plt.legend(loc='upper right', fontsize=10)
	fig.savefig(args.gr_name)
	plt.close(fig)

if __name__ == '__main__':
	args = parse_args()
	L = args.l
	field = init(L) 
	T, Magnetizations, M_theoric = calcul_energy_mag_C_X(field, L)
	plot(args, T, Magnetizations, M_theoric)
	
