import numpy as np

DIR = "/Users/angusr/Python/Subgiants"
n, t2, rms2 = np.genfromtxt("%s/PFS/named_best_times_2_10.txt" % DIR).T
n3, t3, rms3 = np.genfromtxt("%s/PFS/named_best_times_3_10.txt" % DIR).T
n5, t5, rms5 = np.genfromtxt("%s/PFS/named_best_times_5_10.txt" % DIR).T

print len(n)
print n
raw_input('enter')

rms = np.vstack((rms2, rms3, rms5))
t = np.vstack((t2, t3, t5))

dts = []
for i in range(len(n)):
    print n[i]
    print t[0, i], rms[0, i]
    print t[1, i], rms[1, i]
    print t[2, i], rms[2, i]
    l = rms[:, i] == min(rms[:, i])
    print "best rms = ", rms[:, i][l], "best dt = ", t[:, i][l][0]
    print np.where(rms[:, i] == min(rms[:, i]))[0], "\n"
    dts.append(t[:, i][l][0])

np.savetxt("dts.txt", np.array(dts))
