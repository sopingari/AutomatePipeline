
from cc3d.core.PySteppables import *

class WallTest2Steppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        #Code for setting target volumes to initial volumes.
        #Forces cells to roughly maintain constant volume.
        for cell in self.cell_list:
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 1.0
            
            
        self.plot_win1 = self.add_new_plot_window(title='Average Volume',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        self.plot_win2 = self.add_new_plot_window(title='Average Surface Area',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Surface Area', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)  
        self.plot_win3 = self.add_new_plot_window(title='Coziness',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Touching over total', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False) 
        self.plot_win4 = self.add_new_plot_window(title='Compactness',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Compactness', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)  
        self.plot_win1.add_plot("MVol", style='Lines', color='red', size=5)
        self.plot_win2.add_plot("MSur", style='Dots', color='green', size=6)
        self.plot_win3.add_plot("Mcozy", style='lines', color='green', size=4)
        self.plot_win4.add_plot("Comp", style='Lines', color='red', size=5)
 
           
        self.distances = 0
        for CELL1 in self.cell_list_by_type(self.BODY):  #loops over all cells
            #print ("Cell1", "id=", CELL1.id, " type=", CELL1.type)
            for CELL2 in self.cell_list_by_type(self.BODY):  #for each cells, loops over all cells
                #print ("Cell2" "id=", CELL2.id, " type=", CELL2.type)
                vec = self.distance_between_cells(CELL1, CELL2)  #gets the distance to each cell (including itself)
                #print ("vec", vec)
                self.distances += vec #sums the distances
                #print ("looping!  Distance is up to...", self.distances)
        self.mean_distances = self.distances/ ((len(self.cell_list_by_type(self.BODY)))**2)  #once all the distances to each cell is summed, finds the mean
        self.starting_distances = self.mean_distances  # remember this for later
        print ("starting distances")
        print (self.starting_distances)
    

        
            

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        # This if statement allows up to set how often this happens, e.g. every 5 steps
        if mcs % 5 == 0:
            #iterating over all bodies, and summing the volumes and surface areas
            cell_volumes = []
            cell_surfaces = []
            for cell in self.cell_list_by_type(self.BODY):
                cell_volumes.append(cell.volume)
                cell_surfaces.append(cell.surface)
            #dividing the summed volumes and surfaces by the number of cells, to get the average
            mean_cell_volume = sum(cell_volumes)/len(cell_volumes)
            mean_cell_surface = sum(cell_surfaces)/len(cell_surfaces)
        
        
            # Add data to plot arguments are (name of the data series, x, y)
            self.plot_win1.add_data_point("MVol", mcs, mean_cell_volume)
            #print ("id=", cell.id, " type=", cell.type, "volume", cell.volume)
            self.plot_win2.add_data_point("MSur", mcs, mean_cell_surface)      
            #print ("id=", cell.id, " type=", cell.type, "suface", cell.surface)

        if mcs % 5 == 0:
            # "cozy" is the ratio of body-body surface area to total surface area.  
            # It should increase over time due to the favorable body-body contact energy
            cozy = []
            self.touch = 0
            self.notouch = 0
            for cell in self.cell_list_by_type(self.BODY):
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if neighbor.type == self.BODY:
                            self.touch += common_surface_area
                    else: # the neighbor is medium
                        self.notouch += common_surface_area
                coziness = self.touch/(self.touch+self.notouch) #Could probably also make the denominator total surface area; just need to thing about wall effects
                cozy.append (coziness)
            mean_cozy = sum(cozy)/len(cozy)
            # add to the plot
            self.plot_win3.add_data_point("Mcozy", mcs, mean_cozy)
            
            #"comp" (compactness) is a measure of how close the centers of the bodies are to each other relative to when they started
            # it the mean distance of each body to every other body, normalized to the mean distance when the simulation started.  
            # Thus is should start at 1 and generally decrease duing the simulation
            self.distances = 0
            for CELL1 in self.cell_list_by_type(self.BODY):  #loops over all cells
                #print ("Cell1", "id=", CELL1.id, " type=", CELL1.type)
                for CELL2 in self.cell_list_by_type(self.BODY):  #for each cells, loops over all cells
                    #print ("Cell2" "id=", CELL2.id, " type=", CELL2.type)
                    vec = self.distance_between_cells(CELL1, CELL2)  #gets the distance to each cell (including itself)
                    #print ("vec", vec)
                    self.distances += vec #sums the distances
                    #print ("looping!  Distance is up to...", self.distances)
            self.mean_distances = self.distances/ ((len(self.cell_list_by_type(self.BODY)))**2)  #once all the distances to each cell is summed, finds the mean
            compactness = self.mean_distances / self.starting_distances #normalize to how things were to start with
            #print ("compactness", compactness)
            # add to the plot
            self.plot_win4.add_data_point("Comp", mcs, compactness)        

     
    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        # here we specify size of the image saved (1000x1000) - default is 400 x 400
        # resizing of the image is not guaranteed to be implemented

        png_output_path1 = Path(self.output_dir).joinpath("graphOutputVol.png")
        self.plot_win1.save_plot_as_png(png_output_path1, 1000, 1000)
        png_output_path2 = Path(self.output_dir).joinpath("graphOutputSurArea.png")
        self.plot_win2.save_plot_as_png(png_output_path2, 1000, 1000)
        png_output_path3 = Path(self.output_dir).joinpath("graphOutputCozy.png")
        self.plot_win3.save_plot_as_png(png_output_path3, 1000, 1000)
        png_output_path4 = Path(self.output_dir).joinpath("graphOutputCompact.png")
        self.plot_win4.save_plot_as_png(png_output_path4, 1000, 1000)

        # self.plot_win1.save_plot_as_png("uniqueplot1.png", 1000, 1000)
        # self.plot_win2.save_plot_as_png("uniqueplot2.png", 1000, 1000)
        # self.plot_win3.save_plot_as_png("uniqueplot3.png", 1000, 1000)
        # self.plot_win4.save_plot_as_png("uniqueplot4.png", 1000, 1000)
        
        # self.plot_win1.save_plot_as_data("plotdata1.csv")
        
        output_path1 = Path(self.output_dir).joinpath("graphOutputDataVol.txt")
        self.plot_win1.save_plot_as_data(output_path1, CSV_FORMAT)
        output_path2 = Path(self.output_dir).joinpath("graphOutputTestSurArea.txt")
        self.plot_win2.save_plot_as_data(output_path2, CSV_FORMAT)
        output_path3 = Path(self.output_dir).joinpath("graphOutputTestCozy.txt")
        self.plot_win3.save_plot_as_data(output_path3, CSV_FORMAT)
        output_path4 = Path(self.output_dir).joinpath("graphOutputTestCompact.txt")
        self.plot_win4.save_plot_as_data(output_path4, CSV_FORMAT)
        
        
        
        



        