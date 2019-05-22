import os


def submit_0_to_simtime():
    sim_time = {}
    for i in range(1, n_sim):
        print ("____configure file: i", i, "\n")
    
        sim_time[i] = 2*i
    
        f = open(config_file_name, "w")
        f.write(str(0.0)+"\n")
        f.write(str(sim_time[i])+ "\n")
        f.write(str(1)+ "\n")
        f.write(str(0.00002)+ "\n")
        f.write("../../"+ "\n")
        f.write("Data/MatlabGaborFilter/"+ "\n")
        f.write("simpleShapes/"+ "\n")
        f.write("output/synapse_dir/"+ "\n")
        f.close()
    
        os.system("cat "+config_file_name)
    
        os.system("./run_aki_model 1 0 0")
#__
def submit_evey_2_seconds():
    sim_time = 2
    start = 0
    end = 0
    old_synapses_dir = "output/synapse_dir/"
    for i in range(n_sim):
        print ("____configure file: i", i, "\n")
    
        start = end
        end = start + sim_time
        synapses_dir = "output/Start"+str(start)+".000000End"+str(end)+".000000/synapse_dir/"
    
        f = open(config_file_name, "w")
        f.write(str(start)+"\n")
        f.write(str(sim_time)+ "\n")
        f.write(str(1)+ "\n")
        f.write(str(0.00002)+ "\n")
        f.write("../../"+ "\n")
        f.write("Data/MatlabGaborFilter/"+ "\n")
        f.write("simpleShapes/"+ "\n")
        f.write(old_synapses_dir+ "\n")
        f.close()
    
        os.system("cat "+config_file_name)
    
        os.system("./run_aki_model 1 0 0")
    
        old_synapses_dir = synapses_dir

#__

n_sim = 30
config_file_name = "run_config.txt"
submit_evey_2_seconds()
submit_0_to_simtime()
