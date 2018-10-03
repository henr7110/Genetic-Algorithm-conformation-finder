import numpy as np
import timeit
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
R = 8.314*10**(-3) #kJ/molK

def Energy(molecule):
    """returns the sum of energies calculated from angles in the list molecule"""
    return sum(1 + np.cos(w) + np.cos(3*w) for w in molecule)

def Fitness(molecule,T,Z):
    """Returns the boltzmann probability (Fitness) for the molecule at a temp-
    erature T and a partition sum Z"""
    return np.exp(-Energy(molecule)/(R*T)) / Z

def Zeroremover(Fs,k):
    """Takes as input a list of fitness values and an int k specifying the number
    of angles in the molecule. Removes all fitness values that are smaller than
    1/k**3 in order to avoid ZeroDivision exceptions"""
    F = []
    #Itterate through Fitness list
    for i in range(len(Fs)):
        #Only append to returnlist if above 1/N**3
        if Fs[i] > 1.0/k**3:
            F.append(Fs[i])
    return F

def Best(S):
    """Takes as input a list of lists, S, where each element in S is a list of
    angles for that molecule. Returns the molecule from S with lowest energy
    evaluated with the Energy(molecule) function"""
    #Initiate bestEnergy with first element in S
    bestE = Energy(S[0])
    bestmolecule = S[0]
    #Itterate all molecules in S
    for molecule in S[1:]:
        E = Energy(molecule)
        #Check if Energy is smaller
        if E < bestE:
            #If smaller, update best E
            bestE = E
            bestmolecule = molecule
    return bestmolecule, bestE

def Mate(a,b,mutprob):
    """Takes two lists of genotypes: a, b and a float: mutprob. Generates a new
    list: child, with elements chosen randomly from a and b. With a probability
    of mutprob, it then mutates each element in child and returns it."""
    #Generate Child:
    Child = []
    #Itterate through parent elements
    for i in range(len(a)):
        #for each element chose randomly between parent elements and assign
        if np.random.randint(2) == 0:
            Child.append(a[i])
        else:
            Child.append(b[i])
    #Mutate each element in Child:
    for i in range(len(Child)):
        #With probability mutprob, mutate element to new angle
        if np.random.uniform(0,1) < mutprob:
            Child[i] = [np.pi/3, np.pi, np.pi*5/3][np.random.randint(3)]
    return Child

def Simulate(N,K,mutprob,T,stopN,tau=None):
    """Runs genetic algorithm for N molecules of K angles, a mutation probability
    of mutprob, temperature of T, a maximum number of allowed generations stopN
    and a tau if exponential temperature modulation is wanted. Returns the
    number of generations until stopped or converged, as well as the runtime
    and the best molecule configuration for each generation.
    """
    #Initialize variables
    bestEs = []
    Converged = False
    Generation = 0
    t_0 = timeit.default_timer()
    concount = 0

    #generate population:
    S = []
    #Generate N molecules
    for n in range(N):
        molecule = []
        #Append K angles to each molecule
        for i in range(K):
            molecule.append([np.pi/3, np.pi, np.pi*5/3][np.random.randint(3)])
        S.append(molecule)

    #Run genetic algorithm until converged
    while not Converged:
        #Calculate partition sum:
        Z = 0
        for molecule in S:
            Z += np.exp(-Energy(molecule)/(R*T))

        #Calculate Fittness array:
        F = []
        for molecule in S:
            F.append(Fitness(molecule,T,Z))
        #Remove Fitness values below 1/N**3
        F = Zeroremover(F,N)


        #Calculate matingpool
        matingPool = []
        #Find minimum Fitness for matingpool threshold
        mF = min(F)
        #Itterate through fitness-molecule pairs in the population
        for fit,molecule in zip(F,S):
            #Append molecules to the matingpool according to fitness and threshold
            numrep = fit*1/mF
            for n in range(int(numrep)):
                matingPool.append(molecule)

        #mate individuals from matingpool to get newpop
        Newpop = []
        #Generate a new population of N individiuals
        for n in range(N):
            #Randomly select individuals from the matingPool
            a = matingPool[np.random.randint(0,len(matingPool))]
            b = matingPool[np.random.randint(0,len(matingPool))]
            #Append the mated molecules to the new S
            Newpop.append(Mate(a,b,mutprob))
        S = Newpop

        #check for convergence
        #Check if the Generation has reached the maximum allowed
        if Generation == stopN:
            #Save runtime and stop simulation
            dt = timeit.default_timer() - t_0
            Converged = True
            break
        #Check if the best molecule has the optimal energy of -float(K) kJ/mol
        elif Best(S)[1] == -float(K):
            #If so, update convergence counter
            concount += 1
            #If it is the first time, record time of convergence
            if concount == 1:
                dt = timeit.default_timer() - t_0
            #If minima has been reached 10 times in a row, end simulation
            if concount == 10:
                Converged = True
                break
        #If not a minimum, reset counter and keep simulating
        else:
            concount = 0

        #If exponential temperature decay is on, propagate temperature
        if tau != None:
            T = T*np.exp(-1 / tau) + 100

        #Update generation and save lowest energy
        Generation += 1
        bestEs.append(Best(S)[1])
    return Generation, dt , bestEs

def Tempsweeper(sweep,maxG, N, k, mutprob, Tempsweep=True, save=False):
    """Function that sweeps either tau (exponential decay rate of temperature) or
    constant temperature for a given N, k and mutprob. If no Tempsweep or save
    is given, no figure will be saved, and the sweep will be over constant
    temperature, if Tempsweep is set to False, a tau sweep will be done.
    save can be set to the save path for the sweep plot"""
    Gs = []
    #sweep the values presented in list sweep
    for val in sweep:
        #print current value
        print val
        #If tempsweep, input as temperature in simulation
        if Tempsweep:
            Generation, dt,bestEs = Simulate(N, k, mutprob,val,maxG)
        #If tausweep, input as tau in simulation with T_0 = 2400
        else:
            Generation, dt,bestEs = Simulate(N, k, mutprob,2400,maxG,tau=val)
        #Save the sequence of best energies
        Gs.append(bestEs)

    #Plot the data
    fig, ax = plt.subplots(dpi=100)
    #Generate normalized colormap to show different T or Tau
    norm = mpl.colors.Normalize(vmin=sweep[0], vmax=sweep[-1])
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.brg)
    cmap.set_array([])
    #Itterate through BestEs and sweep values (T or tau)
    for E,Ts in zip(Gs,sweep):
        #Plot generations vs. tau/temp and color based on tau/temp
        ax.plot(range(len(E)),E,c=cmap.to_rgba(Ts),linewidth=1)
    #Add colorbar
    cbar = fig.colorbar(cmap, ticks=sweep[::5])
    #Add colorbar label based on temp or tau
    if Tempsweep:
        ax.text(12,-2,'Temperatur (K)',rotation=270)
    else:
        ax.text(10.5,-2.9,'Tau',rotation=270)
    #Set x- and y-labels
    plt.xlabel("Generationer")
    plt.ylabel("Energy (kJ/mol)")
    #If save, save...
    if save != False:
        plt.savefig(save,dpi=700)
    plt.show()

def potentialsurface(surface=True):
    """Function that itterates through two bond angles in the energy-model
    and plots them as either a surface or contour map based on the surface arg-
    ument given in the function definition"""
    #Generate list of angles
    w1_list = np.linspace(0,2*np.pi,300)
    w2_list = np.linspace(0,2*np.pi,300)

    #Generate mesh of energies
    Energysurface = []
    #Itterate through w1
    for w1 in w1_list:
        #For each w1 make a column
        column = []
        for w2 in w2_list:
            #append to column E(w1,w2)
            column.append(Energy([w1,w2]))
        #Build energysurface column by column
        Energysurface.append(column)

    #plot the data based on surface argument
    if surface:
        #Generate omegamesh for 3D plot
        w1s,w2s = np.meshgrid(w1_list,w2_list)
        #generate surface plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(w1s, w2s, np.array(Energysurface), cmap=mpl.cm.coolwarm,
                              linewidth=0, antialiased=False)
        #set axis labels and colorbar with label
        ax.set_xlabel(r"$\omega_1$")
        ax.set_ylabel(r"$\omega_2$")
        colorbar = fig.colorbar(surf, shrink=0.5, aspect=5)
        colorbar.ax.set_ylabel(r"Energy $\frac{kJ}{mol}$")
        plt.show()
    else:
        #plot contour map with angles for x- and y-axis
        plt.imshow(Energysurface,cmap=mpl.cm.coolwarm,extent=[0,2*np.pi,0,2*np.pi])
        #add colorbar with label
        colorbar = plt.colorbar()
        colorbar.ax.set_ylabel(r"Energy $\frac{kJ}{mol}$")
        #label x- and y-axis
        plt.xlabel(r"$\omega_1$")
        plt.ylabel(r"$\omega_2$")
        plt.show()

def boltzmandist():
    configurations = []
    for w_1 in [np.pi/3, np.pi, np.pi*5/3]:
        for w_2 in [np.pi/3, np.pi, np.pi*5/3]:
            for w_3 in [np.pi/3, np.pi, np.pi*5/3]:
                for w_4 in [np.pi/3, np.pi, np.pi*5/3]:
                    for w_5 in [np.pi/3, np.pi, np.pi*5/3]:
                        for w_6 in [np.pi/3, np.pi, np.pi*5/3]:
                            for w_7 in [np.pi/3, np.pi, np.pi*5/3]:
                                configurations.append([w_1,w_2,w_3,w_4,w_5,w_6,w_7])
    temps = range(5,200,75)
    for T in temps:
        #Calculate partition sum:
        Z = 0
        counter = 0
        Energies = [Energy(i) for i in configurations]
        Z = sum(np.exp(-Energy(molecule)/(R*T)) for molecule in configurations)

        p_i = []
        for molecule in configurations:
            p_i.append(Fitness(molecule,T,Z))
        print T
        plt.bar(Energies,p_i,label="T=%dK"%T)
    plt.xlabel(r"$E_i$"+" (kJ/mol)")
    plt.ylabel(r"$p_i$")
    plt.legend()
    plt.savefig("boltzman k=7 T from 5 to 200",dpi=600)
    plt.show()
