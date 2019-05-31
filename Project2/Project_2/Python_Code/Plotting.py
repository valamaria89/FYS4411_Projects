import matplotlib.pyplot as plt


file = open("BF_LR02_SZ1_MC218.txt", 'r')

lines=file.readlines()
iteration = []
energy=[]
for x in lines:
	iteration.append(x.split()[0])
	energy.append(x.split()[1])

iteration = list(map(int, iteration))[1:]
energy = list(map(float, energy))[1:]
#print ([x.split(' ')[1] for x in open(file).readlines()])
#print(iteration)

avg = 0
for i in range(400, len(energy)):
	avg += energy[i]/99

print(avg)
print(len(energy))

plt.plot(iteration, energy)
#plt.ylim([0.4,1])
plt.show()


file.close()