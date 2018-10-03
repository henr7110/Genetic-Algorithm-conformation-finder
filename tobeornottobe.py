import numpy as np
import timeit
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def Genotype(d):
    """Returns the Genotype (list of ascii indexes) from the phenotype (string)
    by conversion of the letters into ascii indexes with the function ord"""
    return(np.array([ord(i) for i in d]))
print Genotype("to be or not to be")

def Phenotype(s):
    """Returns the phenotype (word) from the genotype s which is a list of
    Ascii indexes from 96 to 122 by conversion to string"""
    #Generate word root
    word = ""
    #for Ascii index in s
    for c in s:
        word += chr(c)
    return(word)

def Fitness(a,m):
    """returns normalized fitness from target genotype m given an input genotype
    a by checking how many numbers they have in common"""
    counter = 0
    #itterate through a,m pairs
    for a_i,m_i in zip(a,m):
        #count number of equal elements
        if a_i == m_i:
            counter += 1
    #Return normalized fitness
    return counter / float(len(m))

def Mate(a,b,mutprob):
    """Takes two lists of genotypes: a, b and a float: mutprob. Generates a new
    list: child, with elements chosen randomly from a and b. With a probability
    of mutprob, it then mutates each element in child and returns it."""
    #Generate Child:
    Child = []
    #Itterate through parent elements
    for index in range(len(a)):
        #for each element chose randomly between parent elements and assign
        if np.random.randint(2) == 0:
            Child.append(a[index])
        else:
            Child.append(b[index])

    #Mutate each element in Child:
    for i in range(len(Child)):
        #With probability mutprob, mutate element to new letter
        if np.random.uniform(0,1) < mutprob:
            choice = np.random.randint(96,123)
            #If choice falls on lowest, add a space
            if choice == 96:
                Child[i] = 32
            #Otherwise add the Ascii index directly to the child
            else:
                Child[i] = choice
    return Child

def Simulate(N,mutprob):
    #initialize target string and variables
    mset = np.array(Genotype("to be or not to be"))
    Generation = 1
    t_0 = timeit.default_timer()
    Converged = False

    #generate population
    S = []
    #Generate N sentences
    for i in range(N):
        word = []
        #Append as many letters to each sentence as target sentence
        for d in range(len(mset)):
            choice = np.random.randint(96,123)
            #If choice falls on lowest, add a space
            if choice == 96:
                word.append(32)
            #Otherwise add the Ascii index directly to the sentence
            else:
                word.append(choice)
        S.append(word)

    #Run genetic algorithm until converged
    while not Converged:
        #Get fitness array
        F = []
        for word in S:
            F.append(Fitness(word,mset))
        #Get matingpool
        matingPool = []
        #Itterate through fitness, sentence pairs
        for fit,word in zip(F,S):
            #Append sentences based on fitness
            numrep = fit*100
            for n in range(int(numrep)):
                matingPool.append(word)
        #Mix individuals from mating pool to get newpop
        Newpop = []
        for n in range(N):
            #Randomly select individuals from the matingPool
            a = matingPool[np.random.randint(0,len(matingPool))]
            b = matingPool[np.random.randint(0,len(matingPool))]
            #Append the mated molecules to the new S
            Newpop.append(Mate(a,b,mutprob))
        S = Newpop

        #Check for convergence ie. if the target sentence is in the population
        for word in Newpop:
            #If target sentence reached
            if Fitness(word,mset) == 1.0:
                #Track time, print simulation values and end simulation
                dt = timeit.default_timer() - t_0
                print "Converged in %d Generations and in %.3f " + "seconds with N = %d and mutprob = %.3f" %(Generation,dt,N,mutprob)
                Converged = True
                break
            if Generation == 200:
                dt = timeit.default_timer() - t_0
                print "didn't converge on N= " +str(N) + " mutprob = "+str(mutprob)+ " and dt = %.3f" %(dt)
                Converged = True
                return 200, timeit.default_timer() - t_0
        Generation += 1
    return Generation, dt

def ConvergenceMap(Nlist,mutlist):
    """Function that itterates through all combinations of the values in Nlist
    and mutlist, and returns a textfile with the number of Generations and dt
    for each simulation. Used for constructing a 3D convergencemap of the algorithm"""
    Glist = []
    Tlist = []
    #Itterate through Nlist
    for N in Nlist:
        #Itterate through mutlist
        for mutprob in mutlist:
            #run simulation with parameters
            Generations, dt = Simulate(N,mutprob)
            Glist.append(Generations)
            Tlist.append(dt)
    string = "N,mutprob,Generations,Time\n"
    #Generate output string with parameters
    for N,n1 in zip(Nlist,range(len(Nlist))):
        for mutprob,n2 in zip(mutlist,range(len(mutlist))):
            string += str(N)+","+str(mutprob)+","+
                      str(Glist[(n1*len(mutlist))+n2]) + "," +
                      str(Tlist[(n1*len(mutlist))+n2])+"\n"

def Convergencereader(path):
    """Function that reads the print output from ConvergenceMap and plots the 3D
    Convergence surface in time an generations """
    file = open(path,"r")
    N = []
    mutprob = []
    Generations = []
    time = []
    for line in file:
        if line[0] == "d":
            #did't converge
            s = line.split("= ")
            mutprob.append(float(s[2]))
            N.append(int(s[1].split(" ")[0]))
            Generations.append(200)
            time.append(200.0)
        else:
            #did converge
            b = line.split(" ")
            Generations.append(float(b[2]))
            time.append(float(b[5]))
            N.append(int(b[15]))
            mutprob.append(float(b[-1]))

    mp = np.reshape(mutprob,(12,50)).T
    n = np.reshape(N,(12,50)).T
    g = np.reshape(Generations,(12,50)).T
    t = np.reshape(time,(12,50)).T
    max(time)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(mp, n, t, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_xlabel("mutprob")
    ax.set_ylabel("N")
    ax.set_zlabel("Time until converged")
    # plt.savefig("Convergencemap",dpi=500)
    plt.show()

def allinone(Ns,string,k,save=False):
    """Function that generates a population and checks whether all letters in
        a target sentence (string) are represented somewhere in the population.
        Only k possible letters (from a to space) are used where k can range from
        1 to 23. Afterwards it plots a histogram over the fraction of target
        letters contained in the population based on population length.
        These N-sweeps are specified from the Ns list. Can be saved by passing
        a path to the save argument."""
    in0 = []
    #For each population size in Ns
    for N in Ns:
        #print the population size
        print N
        #Initialze counters
        allinpop = 0
        notinpop = 0
        #Run 500 populations of size N
        for i in range(500):
            #Print runnumber for every 100
            if i%100==0:
                print i
            #Convert string to Ascii indexes
            mset = Genotype(string)

            #generate population
            S = []
            #Generate N sentences
            for i in range(N):
                word = []
                #Append as many letters to each sentence as target sentence
                for d in range(len(mset)):
                    choice = np.random.randint(96,96 + k)
                    #If choice falls on lowest, add a space
                    if choice == 96:
                        word.append(32)
                    #Otherwise add the Ascii index directly to the sentence
                    else:
                        word.append(choice)
                S.append(word)

            #Check if all target letters are represented at correct positions in
            #the population
            correct = np.zeros(len(mset))
            #Itterate through all sentences in S
            for sentence in S:
                #Subtract target sequence from sentence
                diff = sentence-mset
                #Check diff values for zeroes
                for letter,num in zip(diff,range(len(diff))):
                    #If the subtraction yields 0 somewhere, that was a correct letter
                    #If this isn't already in the correct list, change it to 1
                    if letter == 0 and correct[num] == 0:
                        correct[num] = 1
            #If all indicies are represented in the population up allinpop,
            #otherwise up notinpop
            if sum(correct) == len(mset):
                allinpop += 1
            else:
                notinpop += 1
        #Save data
        in0.append(float(allinpop)/(notinpop+allinpop))

    #Plot results
    p1 = plt.bar(Ns, in0, width=1,)
    p2 = plt.bar(Ns, 1-np.array(in0),width=1, bottom=in0)
    plt.ylabel("Andel")
    plt.xlabel("Startpopuplation")
    plt.xticks(Ns)
    plt.legend((p1[0],p2[0]),("Alle korrekt","Ikke alle korrekt"))
    #Save plot
    if save != False:
        plt.savefig(save,dpi = 600)
    plt.show()
