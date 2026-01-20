#define _USE_MATH_DEFINES

#define MAX_NEIGHBOURS    4
#define MAX_SPECIES       150
#define MAX_CELLS         1000

#define MAX_REACTIONS     500
#define MAX_REACTION_SIZE 4

const int MAX_EVENTS = MAX_CELLS*(MAX_REACTIONS + MAX_NEIGHBOURS*MAX_SPECIES);



#include <cstdlib>
#include <cmath>
#include <iostream>     // For basic input output
#include <iomanip>      // IO manipulation, setting width
#include <fstream>      // For opening files

#include <string>
#include <sstream>      // For parsing model text files



#include <random>

// #include <SFML/Window.hpp>
// #include <SFML/Graphics.hpp>
// #include <SFML/System.hpp>

using namespace std;

void tell(string message, bool fail){
    cout << message << endl;
    cin.get();
    if(fail){
        exit(1);
    }
}


// Restructuring:
// We make a Reaction class with the species index, state change, number of reactants, reaction name etc.
// Hence we have a Reaction for each actual reaction, and one corresponding to the diffusion of each species.
// We have an Event class containing only target cells and an index or pointer referencing the original Reaction that it refers to.



// Generic reaction
class Reaction {

    public:

        // TODO rename species to reaction_species to avoid confusion with Model::species.
        // Fixed properties of the reaction
        string reaction_name;
        int species[MAX_REACTION_SIZE]; //Declares and array of [] integers, MAX_REACTION_SIZE, in this case.
        int state_change[MAX_REACTION_SIZE];

        // Total number of nonunique species either created or destroyed. I.e. number of state changes that must be executed
        int n_reaction_species;

        // The species indices that react and determine the propensity. Index is the index in species list in reactants.
        int reactants[MAX_REACTION_SIZE];
        int n_reactants;

        // Rate. In case reaction is a diffusion, it is the diffusion coefficient
        double k;

        // Number of events of this reaction that have occured
        int count;

        // Constructor
        Reaction();
};

Reaction::Reaction() {

    n_reaction_species = 0;
    n_reactants = 0;

    count = 0;
}

class Event {

    public:
        int target_cells[MAX_REACTION_SIZE];
        double rate;
        Reaction* reaction_ptr;

        // 0 is diffusion, 1 is reaction.
        
};



class Cell {

    public:

        // ID
        int idx;

        //time
        //double time;

        // Positions
        double x;
        double y;

        // Neighbours
        int neighbours[MAX_NEIGHBOURS];
        double neighbour_dist[MAX_NEIGHBOURS];
        int n_neighbours;
        
        // Current number of each species in the cell
        double state[MAX_SPECIES];


        // Volume and area for correcting the propensities
        double volume;
        double area;

        Cell();
        void set_species_number(int species_idx, long long int value);
        void add_neighbour(int neighbour_idx, double distance);

};

Cell::Cell(){

    n_neighbours = 0;

    for(int i = 0; i < MAX_SPECIES; i++){
        state[i] = 0;
    }
}


void Cell::set_species_number(int species_idx, long long int value){

    // Add check that species exists

    state[species_idx] = int (value);
}

void Cell::add_neighbour(int neighbour_idx, double distance){

    // Add overflow check
    if(n_neighbours == MAX_NEIGHBOURS){
        cout << "Number of neighbours of cell " << idx << " exceeds maximum number of neighbours" << endl;
        cin.get();
        exit(1);
    }
    neighbours[n_neighbours] = neighbour_idx;
    neighbour_dist[n_neighbours] = distance;
    n_neighbours++;

}

// Random number generation
random_device rd;  							// Will be used to obtain a seed for the random number engine. As rd() not guaranteed to follow any distribution.
mt19937 gen(rd());							// Standard mersenne_twister_engine seeded with rd(). Is statistically random.
uniform_real_distribution<> dis(0.0, 1.0); 	// Call to "dis(gen)" will now generate a random double in [0,1)


class Model {

    public:

        // Name of the model
        string model_name;

        // Bool to check if initial recording has already happened
        bool initial_recording;

        // Time and iteration
        double t;
        long long int i_iter;

        // The cells that make up the grid
        // TODO: Rename grid to cells or cell_list. Is not always a grid.
        Cell grid[MAX_CELLS];
        int n_cells;

        // Variables to allow rendering
        double min_x, max_x, min_y, max_y, min_dist;

        // The possible reactions that can occur
        Reaction reactions[MAX_REACTIONS];
        int n_reactions;

        Reaction diffusion_reactions[MAX_SPECIES];

        // The chemical species
        string species[MAX_SPECIES];
        int species_dimensions[MAX_SPECIES];
        int n_species;

        // All possible reactions and their propensities
        Event events[MAX_EVENTS];
        int event_count[MAX_REACTIONS + MAX_SPECIES];
        int n_events;
        double propensities[MAX_EVENTS];
        double a0;

        Reaction* reaction_ptr; // Reaction pointer for temporary access to reactions.
        Event* event_ptr;       // Event pointer for temporary access to events.
        Cell* cell_ptr;         // Cell pointer for temporary access to cells.

        
        // Constructor
        Model();

        // TODO: create new function which assigns reactants. Now there is code duplication in add_species and add_reaction.
        // TODO: add datestring to logging file names. Prevents accidentally overwriting data.

        void create_cartesian_grid(double d_x, double d_y, double d_z, int n_x, int n_y, bool periodic);

        void run(long long int iterations, long long int break_iter, bool verbose, string folder_name);
        void load_model(string file_name);
        void load_initial_conditions(string file_name);
        // void set_receptors(int low_power, int high_power, int n_cells);


        void record_concentrations(string project_name, double n_iter);

        void record_all_species_in_cell(string project_name, const double n_iter, const int i_cell);
        void record_all_cells_for_species(string project_name, string species_name,int* cur_cells, const double n_iter);//function to get graphs


    private:

        void add_species(string species_name, int dimension, double dif_coefficient);
        void add_reaction(string reaction_species[], int state_changes[], double rate, int n_reaction_species, string reaction_name);

        int get_species_idx(string species_name);

        void calculate_propensities();//int i_cell);
        void update_state(int mu, double dt);//, int i_cell);

};

Model::Model(){
    n_cells = 0;
    n_reactions = 0;
    n_species = 0;
    n_events = 0;

    t = 0;

    min_x = min_y = 0;
    max_x = max_y = 1;

    initial_recording = true;

    for(int i_event = 0; i_event < MAX_EVENTS; i_event++){
        propensities[i_event] = 0;
    }
}


void Model::record_all_species_in_cell(string project_name, const double n_iter, const int i_cell) {
	
    // Recording the number of components as a function of time in a given cell.
    
    // Setting file name to include project and cell index.
    string file_name, cell_name;
    cell_name = "_components_in_cell_" + to_string(i_cell);
    file_name = project_name + cell_name + ".dat";

	ofstream datafile;

	if(initial_recording) {
		
		datafile.open(file_name);
		//datafile << "#Evolution of concentrations of " << model_name << " species in one grid cell" << endl;
		datafile << "#Evolution of concentrations of the species in one grid cell" << endl;
		datafile << "#The components on grid cell " << i_cell << endl << endl;
		datafile << setw (20) << "# Step (n)" << setw (20) << "Time (t)";

        for(int i_species = 0; i_species < n_species; i_species++){
            datafile << setw(20) << species[i_species];
            // datafile << "," << species[i_species];
        }

		datafile << endl << endl;
				
	}
	
	datafile.open(file_name, ios::app);
	
	datafile << setw (20) << n_iter << setw (20) << t; //cell_ptr->time;
	for (int i_species = 0; i_species < n_species; i_species++) {
		datafile << setw(20) << grid[i_cell].state[i_species];
	}

	datafile << endl;
	datafile.close();
}

void Model::record_all_cells_for_species(string project_name, string species_name, int* cur_cells, const double n_iter) {
	
    // Recording the number of components as a function of time in a given cell.
    
    // Setting file name to include project and cell index.
    string file_name, species_description;
    species_description = "_concentration_of_species_" + species_name;
    file_name = project_name + species_description + ".dat";

	ofstream datafile;

	if(initial_recording) {
		
		datafile.open(file_name);
		datafile << "#Evolution of concentrations of " << species_name << " in " << model_name << " in all cells" << endl;
        cout << endl;

		datafile << "\n\n\n" << setw (20) << "# Step (n)" << setw (20) << "Time (t)";

        for(int i_cell = 0; i_cell < n_cells; i_cell++){
            datafile << setw(18) << "Cell " << setw(2) << i_cell;
        }
		datafile << endl << endl;
				
	}
	
	datafile.open(file_name, ios::app);
	
	datafile << setw (20) << n_iter << setw (20) << t;
	for (int i_cell = 0; i_cell < n_cells; i_cell++) {
		datafile << setw(20) << cur_cells[i_cell];
	}

	datafile << endl;
	datafile.close();
}

void Model::record_concentrations(string project_name, double n_iter) {

    int species_distribution[MAX_CELLS];

    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        record_all_species_in_cell(project_name, n_iter, i_cell);
    }

    for(int i_species = 0; i_species < n_species; i_species++){
        for(int i_cell = 0; i_cell < n_cells; i_cell++){
            species_distribution[i_cell] = grid[i_cell].state[i_species];
        }

        record_all_cells_for_species(project_name, species[i_species], species_distribution, n_iter);
    }

    initial_recording = false;

}

void Model::load_initial_conditions(string file){

    // TODO IDEA. Allow marking of cells with flags. Then add species randomly or uniformly over those flags.

    // TODO: code duplication, make a function
    string folder = "";
    string path;
    path = folder + file;

    ifstream setup_file;
    setup_file.open(path);
    if (!setup_file){
        cout << "File " << path << " does not exist" << endl;
        cin.get();
        exit(1);
    }

    string line, species_name;
    stringstream stream;

    // Loaded from file
    int i_species, count;

    // For random distribution of particles.
    double p_cell, r;
    int i_cell;

    if(n_cells == 0){
        tell("Number of cells must be greater than 0", 1);
    }

    p_cell = 1/double(n_cells);


    cout << "Initial conditions: " << endl;
    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        stream.clear();
        stream << line;
        stream >> species_name >> count;

        cout << species_name << ": " << count << endl; //count = number of particles in this case

        i_species = get_species_idx(species_name);

        // for(int i_particle = 0; i_particle < count; i_particle++){
        //     r = dis(gen);
        //     i_cell = int(r/p_cell);
        //     // cout << "Cell: " << i_cell << endl;
        //     grid[i_cell].state[i_species] += 1;
        // }
        for (int i_cell = 0; i_cell<n_cells; i_cell++){
            grid[i_cell].state[i_species] = count/n_cells;
        }
    }

}

void Model::load_model(string file){
    
    string folder = "";
    string path;
    path = folder + file;

    

    ifstream setup_file;
    setup_file.open(path);
    if (!setup_file){
        cout << "File " << path << " does not exist" << endl;
        cin.get();
        exit(1);
    }


    // Variables for data parsing
    stringstream stream;
    string line;

    string species_name, reaction_name;
    float rate;
    int species_dimensions;
    int indv_change;
    int n_reaction_species;
   
    string reactants[MAX_REACTION_SIZE];
    int state_change[MAX_REACTION_SIZE];


    // Loading in model name
    while (!setup_file.eof()) {

        getline(setup_file, line);

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        model_name += line;
    }

    // Loading in species data
    while (!setup_file.eof()) {

        // Save the line in "line".
        getline(setup_file, line); 

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // If symbol indicates next section, break this loop
        if(line[0] == '>'){
            break;
        }

        // Tell stream that it has not reached the end.
        stream.clear();

        stream << line;
        stream >> species_name >> species_dimensions >> rate;
        cout << "Species " << species_name << ", dimension: " << species_dimensions << ", rate: " << rate << endl;
        cout << endl;

        // Sanitizing species name:
        string forbidden_chars = "<>:\"/\\|?*";
        for(int i_char = 0; i_char < forbidden_chars.length(); i_char++){
            if (species_name.find(forbidden_chars[i_char]) != string::npos) {
                cout << "Species name contains forbidden character: " << forbidden_chars[i_char] << endl;
                cin.get();
                exit(1);
            }
        }

        add_species(species_name, species_dimensions, rate);
    }

    // Loading in reaction data
    while (!setup_file.eof()) {

        // Save the line in "line".
        getline(setup_file, line); 

        // If comment or empty line, skip line
        if(line[0] == '/' || line.length() == 0){
            continue;
        }

        // Get reaction name
        reaction_name = line;
        cout << reaction_name << endl;

        getline(setup_file, line); 

        stream.clear();

        // Get reactants
        stream << line;

        n_reaction_species = 0;
        while (!stream.eof()){
            stream >> species_name >> indv_change;

            cout << "Species " << species_name << " , change: " << indv_change << endl;
            reactants[n_reaction_species] = species_name;
            state_change[n_reaction_species] = indv_change;
            n_reaction_species++;
        }
        // n_reaction_species--; //BOQI CHANGED, was not commented out before. But commenting it out fixed the error of a receptor and unbound particle dissapearing without a bound particle appearing (leading to depletion of receptors)

        // Get reaction rate
        stream.clear();
        getline(setup_file, line);
        stream << line;
        stream >> rate;
        cout << "Rate " << rate << endl;
        cout << endl;

        add_reaction(reactants, state_change, rate, n_reaction_species, reaction_name);

        // TODO: add check that some combination of reactants remain constant.
    }
}


// Mostly done
void Model::add_species(string species_name, int dimension, double dif_coefficient){

    if(n_species == MAX_SPECIES){
        tell("Number of species exceeds maximum number of species", 1);
    }

    // Safety check: diffusion coefficient cannot be negative
    if(dif_coefficient < 0){
        cout << "Diffusion coefficient for species " << species_name << " must be greater or equal to 0";
        cin.get(); exit(1);
    }

    // Check if species exists already
    for(int i = 0; i < n_species; i++){
        if(species_name == species[i]){
            cout << "Species already exists";
            return;
        }
    }

    species[n_species] = species_name;
    species_dimensions[n_species] = dimension;

    reaction_ptr = &diffusion_reactions[n_species];
    reaction_ptr->k = dif_coefficient;

    // In both cells, the same species is involved
    reaction_ptr->species[0] = n_species;
    reaction_ptr->species[1] = n_species;
    reaction_ptr->n_reaction_species = 2;

    // In one cell, a species disappears, which appears in the other cell
    reaction_ptr->state_change[0] = -1;
    reaction_ptr->state_change[1] = 1;

    // Add reactants
    reaction_ptr->n_reactants = 1;
    reaction_ptr->reactants[0] = 0;

    reaction_ptr->reaction_name = species_name + " diffusion";

    n_species++;
}

int Model::get_species_idx(string species_name){

    int species_idx = -1;

    for(int i = 0; i < n_species; i++){
            if(species[i] == species_name){
                species_idx = i;
            }
    }

    if(species_idx == -1){
        cout << "Reaction species " << species_name << " was not yet defined" << endl;
        cin.get();
        exit(1);
    }

    return species_idx;
}

void Model::add_reaction(string reaction_species[], int state_changes[], double rate, int n_reaction_species, string reaction_name){

    if(n_reactions == MAX_REACTIONS){
        tell("Number of reactions exceeds maximum number of reactions", 1);
    }

    if(n_reaction_species > MAX_REACTION_SIZE){
        cout << "Number of reactants in reaction " << n_reactions + 1 << " exceeds maximum number of reactants" << endl;
        cin.get();
        exit(1);
    }

    // Safety check: reaction rate must be greater or equal to 0
    if(rate < 0){
        cout << "Reaction rate coefficient must be greater or equal to 0";
        cin.get(); exit(1);
    }

    // Get pointer to current reaction
    reaction_ptr = &reactions[n_reactions];
    int species_idx;
    
    // Loop over all species in the reaction to see if they exist, and if so, add their index to the reaction class
    for(int i = 0; i<n_reaction_species; i++){
        
        species_idx = get_species_idx(reaction_species[i]);

        reaction_ptr->species[i] = species_idx;
        reaction_ptr->state_change[i] = state_changes[i];

        // Determining which species are reactants
        if(state_changes[i] < 0){
            reaction_ptr->reactants[reaction_ptr->n_reactants] = i;
            reaction_ptr->n_reactants++;
        }
        
    }

    reaction_ptr->k = rate;
    reaction_ptr->n_reaction_species = n_reaction_species;

    reaction_ptr->reaction_name = reaction_name;
    
    n_reactions++;

}

// Mostly done
void Model::create_cartesian_grid(double d_x, double d_y, double d_z, int n_x, int n_y, bool periodic){

    if(n_x*n_y > MAX_CELLS){
        tell("Number of cells exceeds maximum number of cells", 1);
    }
    // if(n_x <= 1 || n_y <= 1){
    //     tell("Grid must have a 2 or more points in both dimensions", 1);
    // }

    if(d_x <= 0 || d_y <= 0 || d_z <= 0){
        tell("All lengths must be nonzero and positive", 1);
    }

    double dist_x, dist_y;
    dist_x = d_x/(n_x-1);
    dist_y = d_y/(n_y-1);

    max_x = d_x;
    max_y = d_y;
    min_dist = min(dist_x, dist_y);

    n_cells = n_x*n_y;

    // Creating periodic boundary conditions
    for(int i = 0; i < n_cells; i++){
        grid[i].idx = i;
        grid[i].x = double(i%n_x)*dist_x;
        grid[i].y = double(i/n_x)*dist_y;

        grid[i].area   = (d_x/n_x)*(d_y/n_y);
        grid[i].volume = grid[i].area * d_z;

        if(i%n_x != 0) {        // Cells not at left boundary
			grid[i].add_neighbour(i - 1, dist_x);
		} else if(periodic) {   // Cells at left boundary
			grid[i].add_neighbour(i + n_x - 1, dist_x);
		}
		
		if(i%n_x != n_x - 1) {  // Cells not at right boundary
			grid[i].add_neighbour(i + 1, dist_x);
		} else if(periodic) {   // Cells at right boundary
			grid[i].add_neighbour(i - n_x + 1, dist_x);
		}
		
		if(i > n_x - 1) {       // Cells not at top boundary
			grid[i].add_neighbour(i - n_x, dist_y);
		} else if(periodic) {   // Cells at top boundary
			grid[i].add_neighbour(i + n_cells - n_x, dist_y);
		}
		
		if(i < n_cells - n_x) { // Cells not at bottom boundary
			grid[i].add_neighbour(i + n_x, dist_y);
		} else if(periodic) {   // Cells at bottom boundary
			grid[i].add_neighbour(i%n_x, dist_y);
		}
    }
}

void Model::run(long long int iterations, long long int break_iter, bool verbose, string folder_name){

    // Temporary variables.
    Cell cell;
    int i_species, species_dim;
    
    //! BOQI added
    string file_name = folder_name + string("parameters") + ".dat";

	ofstream datafile;
    datafile.open(file_name, ios::app);
    datafile << "Number of arms" << setw(20) << n_species-7 << endl;
    datafile.close();
    //! 

    

    // Define events
    
    // Diffusion
    
    for(int i_cell = 0; i_cell < n_cells; i_cell++){
        cell = grid[i_cell];

        for(int i_neighbour = 0; i_neighbour < cell.n_neighbours; i_neighbour++){

            for(int i_species = 0; i_species < n_species; i_species++){

                // New code
                event_ptr = &events[n_events];
                reaction_ptr = &diffusion_reactions[i_species];

                event_ptr->reaction_ptr = reaction_ptr;
                event_ptr->target_cells[0] = i_cell;
                event_ptr->target_cells[1] = cell.neighbours[i_neighbour];
                event_ptr->rate = (reaction_ptr->k)/pow(cell.neighbour_dist[i_neighbour], 2);

                n_events++;
                // End new code
            }
        }
    }

    // Reactions
    
    cout << "Number of reactions: " << n_reactions << endl;
    cout << "Number of cells: " << n_cells << endl;
    

    cout << endl;
    cout << endl;

    for(int i_cell = 0; i_cell < n_cells; i_cell++){

        for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){

            // Get data from reaction class for species and rate
            reaction_ptr = &reactions[i_reaction];   

            event_ptr = &events[n_events];
            event_ptr->reaction_ptr = reaction_ptr;

            for(int i_reaction_species = 0; i_reaction_species < reaction_ptr->n_reaction_species; i_reaction_species++){
                event_ptr->target_cells[i_reaction_species] = i_cell;
            }

            // TODO add in more documentation for rate correction.
            // When adding reaction, need to add which right hand reagent k is relative to.
            // In file, denote this by a new line.
            
            // Set base event rate
            event_ptr->rate = reaction_ptr->k;

            // Correct for volume and area
            // All constants are relative to the membrane species.
            event_ptr->rate *= grid[i_cell].area;

            cout<<"area is "<<grid[i_cell].area<<endl;
            cout<<"volume is "<<grid[i_cell].volume<<endl;

            for(int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){
                i_species = reaction_ptr->species[reaction_ptr->reactants[i_reactant]];
                species_dim = species_dimensions[i_species];

                if(species_dim == 2) {
                    event_ptr->rate /= grid[i_cell].area;
                } else if(species_dim == 3) {
                    event_ptr->rate /= grid[i_cell].volume;
                } else {
                    cout << "The species dimension " << species_dim << "is not 2 or 3. Check the input file." << endl;
                    cin.get();
                    exit(1);
                }
            }
            
            // Update event count
            n_events++;

        }
    }

    double a1 = 0;
    double r1, r2, dt;
    int mu;

    // Eventually replace by while loop!!!
    for(i_iter = 0; i_iter < iterations; i_iter++){

	if (i_iter%1000000==0) cout << "i_iter = " << i_iter << endl;
        // // Cdc42 test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += grid[i_cell].state[0];
        //     const_count += grid[i_cell].state[3];
        //     const_count += grid[i_cell].state[4];
        //     const_count += grid[i_cell].state[6];

        // }
        // cout << "Total cdc42: " << const_count << endl;

        // if (const_count != 1000){
        //     cout << "Event: " << events[mu].reaction_ptr->reaction_name << endl;
        //     cin.get();
        // }

        // // Bem1 test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += grid[i_cell].state[1];
        //     const_count += grid[i_cell].state[7];
        //     const_count += grid[i_cell].state[8];

        // }
        // cout << "Total Bem1: " << const_count << endl;

        // GEF test.
        // const_count = 0;
        // for(int i_cell = 0; i_cell < n_cells; i_cell++){
        //     const_count += grid[i_cell].state[2];
        //     const_count += grid[i_cell].state[8];
        //     const_count += grid[i_cell].state[9];

        // }
        // cout << "Total GEF: " << const_count << endl;

        // Recalculate propensities.
        calculate_propensities();

        // for(int i_cell = 0; i_cell < n_cells; i_cell++){

        //     // Recalculate propensities.
        //     calculate_propensities(i_cell);

        //     // Picking time of reaction	
        //     // Evaluates the time before the next reaction will take place, storing it in dt (Eq. 21a in Gillespie 1977)
        //     cell_ptr = &(grid[i_cell]);
        //     r1 = dis(gen); //random number in homogeneous distribution
        //     dt = -log(r1)/a0;
        //     cell_ptr->time += dt;

        //     // Picking event
        //     // Picks the reaction that is going to take place (Eq. 21b in Gillespie 1977)
        
        //     a1 = 0;
        //     r2 = dis(gen); //random number in homogeneous distribution
            
        //     for (int i_reaction = 0; i_reaction < n_reactions; i_reaction++) {

        //         a1 += propensities[i_reaction];
        //         if (r2*a0 <= a1){
        //             mu = i_reaction;
        //             break;
        //         }
        //     }

        //     // Updating the state vector with the chosen reaction
        //     update_state(mu, dt, i_cell);
        // }

        //original
        r1 = dis(gen);
        dt = -log(r1)/a0;
        t += dt; //BOQI this is where the time difference is calculated
        

        // Picking event
        // Picks the reaction that is going to take place (Eq. 21b in Gillespie 1977)
        
        a1 = 0;
        r2 = dis(gen);
        
        for (int i_event = 0; i_event < n_events; i_event++) {

            a1 += propensities[i_event];
            if (r2*a0 <= a1){
                mu = i_event;
                break;
            }
        }
	
        
        if(i_iter%break_iter == 0){ // and i_iter/break_iter>9*iterations/(10*break_iter)

            if(verbose){

                cout << "Iteration number " << i_iter << endl;

                for(int i_event = 0; i_event < n_events; i_event++){
                    
                    reaction_ptr = events[i_event].reaction_ptr;
                    cout << i_event << " " << reaction_ptr->reaction_name << "Cell " << events[i_event].target_cells[0] << ": " << propensities[i_event] << "   " << endl;
                }
                cout << "Mu = "<< mu << endl;

                for(int i_cell = 0; i_cell < n_cells; i_cell++){
                    for(int i_species = 0; i_species < n_species; i_species++){
                        cout << grid[i_cell].state[i_species] << "   ";
                    }
                    cout << endl;
                }
                cin.get();
                
            }
            //This is where the output is stored
            record_concentrations(folder_name + string("Testrun_" + to_string(n_species-7) + string("_arms")), double(i_iter));
            //BOQI changed
        }


        // Updating the state vector with the chosen reaction
        update_state(mu, dt);

    }
}

void Model::calculate_propensities(){//int i_cell){

    double propensity;
    int i_reaction_species, i_species;
    int reactant_number;

    // Reset total propensity
    a0 = 0;

    // Loop over all events, and calculate propensities for each event.
    for(int i_event = 0; i_event < n_events; i_event++){

        event_ptr = &events[i_event];
        reaction_ptr = event_ptr->reaction_ptr;
        propensity = event_ptr->rate;

    //     // Loop over all events, and calculate propensities for each event. changed events to reaction to only do it for one cell
    // for(int i_reaction = 0; i_reaction < n_reactions; i_reaction++){

    //     reaction_ptr = &reactions[i_reaction];
    //     // reaction_ptr = event_ptr->reaction_ptr;
    //     propensity = reaction_ptr->k;//k=rate of reaction

        // cout << i_event << " " << reaction_ptr->reaction_name << " in cell " << events[i_event].target_cells[0] << endl;
        // cout << "Involved reactants and their number: " << endl;
        
        for(int i_reactant = 0; i_reactant < reaction_ptr->n_reactants; i_reactant++){

            // Index "i_reactant" is an index over all reactants in the reaction.
            // It points to an entry in "reactants" which contains an index i_reaction_species.
            // This points to an entry in "species" in Model which is the index of the actual reaction.
            i_reaction_species = reaction_ptr->reactants[i_reactant];
            i_species = reaction_ptr->species[i_reaction_species];
            
            // Look at target cell, and in that cell, get number of molecules of the species.
            reactant_number = grid[event_ptr->target_cells[i_reaction_species]].state[i_species];
            // reactant_number = grid[i_cell].state[i_species];//search
            propensity *= double(reactant_number);

            // cout << i_species << " " << species[i_species] << ": ";
            // cout << reactant_number << endl;
            // cin.get();

        }


        

        // Safety check: each individual propensity must be positive
        if(propensity < 0){
            cout << "Propensity for event " << i_event << " went negative";
            cin.get(); exit(1);
        }

        propensities[i_event] = propensity;
        a0 += propensity;
        
    }
}

void Model::update_state(int mu, double dt){//, int i_cell){


    event_ptr = &events[mu];
    reaction_ptr = event_ptr->reaction_ptr;
    // reaction_ptr = &reactions[mu];
    // reaction_ptr = event_ptr->reaction_ptr;

    // Count the reaction
    reaction_ptr->count++;

    for(int i_reactant = 0; i_reactant < reaction_ptr->n_reaction_species; i_reactant++){
        
        // Get pointer to target cell
        cell_ptr = &(grid[event_ptr->target_cells[i_reactant]]);
        // cell_ptr = &(grid[i_cell]);

        // Update state of target cell
        cell_ptr->state[reaction_ptr->species[i_reactant]] += reaction_ptr->state_change[i_reactant];

        // Safety check: population of state must always be 0 or bigger
        if(cell_ptr->state[reaction_ptr->species[i_reactant]] < 0){
            cout << "Population of cell  " << cell_ptr->idx << ", reactant " << i_reactant << " went negative";
            cin.get(); exit(1);
        }
    }
}

vector<vector<double>> create_colour_scale(int n_colours, int max_n_species){

    vector<vector<double>> colour_scales(n_colours, vector<double>(max_n_species, 0));

    double x_mid, width;
    double x;

    for(int i_colour = 0; i_colour < n_colours; i_colour++){

        x_mid = i_colour * 1.0/(n_colours-1);
        width = 1.0/(n_colours-1);


        for(int i = 0; i < max_n_species+1; i++){

            // x = i/double(max_n_species);
            x = log(i+1)/log(max_n_species+1);
            
            if ((x < x_mid - width) || (x > x_mid + width)){
                colour_scales[i_colour][i] = 0;
            } else {
                colour_scales[i_colour][i] = 0.5 * (cos(M_PI  * (x-x_mid)/(width)) + 1);
            }
        }
    } 

    return colour_scales;
}

// void render(Model* model_ptr){

//     // In creating grid, store min and max x and y values. Use these to 
//     // determine x scale and y scale, use smallest scale to draw.

//     // Internal coordinates are in pixels. Starting at top left.
//     // Attachment point of shapes (rectangle, circle) is at top left.

//     Reaction* reaction_ptr;
//     Cell* cell_ptr;

//     double x0, y0;
//     int plot_species = 0;
//     string LocalComposition;

//     // For time printing
//     sf::Clock clock;
//     sf::Time Current;
//     // For cycling through species
//     sf::Time LeftClick = clock.getElapsedTime(), RightClick = clock.getElapsedTime();  


//     // Determining size of drawn images
//     int width = 1600;
//     int height = 900;

//     sf::RenderWindow window(sf::VideoMode(width, height), "Reaction-Diffusion");

//     // Area of drawable canvas does not include the boder
//     double border = 50.f, s_border = 2.5;
    
//     // For printing text;
//     int text_size = 18;
//     double line_size = 20;

//     // Idea: the cells in the model have certain coordinates. We use these coordinates to draw a map of the network. 
//     // We determine the size scaling by comparing the maximum distance in internal coordinates in the x and y direction,
//     // with the max size of the window (minus a border size).
    
//     double v_scale, h_scale, scale, s_size, l_size;
//     width -= 2*border;
//     height -= 2*border;

//     // Total distance in internal coordinates
//     double tot_x, tot_y;
//     tot_x = (model_ptr->max_x - model_ptr->min_x) + model_ptr->min_dist;
//     tot_y = (model_ptr->max_y - model_ptr->min_y) + model_ptr->min_dist;

//     v_scale = width/tot_x;
//     h_scale = height/tot_y;
//     scale = min(v_scale, h_scale);

//     // Debug statements
//     // cout << "Max x: " << model_ptr->max_x << endl;
//     // cout << "Min dist: " << model_ptr->min_dist << endl;  
//     // cout << "Totx: " << tot_x << endl;
//     // cout << "Toty: " << tot_y << endl;
//     // cout << "Scale: " << scale << endl;
//     // cin.get();

//     // The maximum size individual cells can be drawn at, without them touching.
//     s_size = scale*(model_ptr->min_dist) - 2*s_border;


//     // Legenda variables.
//     double count_species, factor;
//     double colorscale_R, colorscale_G, colorscale_B;

//     int n_colours = 3;
//     int N_legend = 15, max_n_species = 10e7;//check high_power if the power is the same

//     // The size of the legenda boxes.
//     l_size = double(height)/(N_legend) - 2*s_border;
    

//     // Color scale
//     vector<vector<double>> colour_scale;
//     colour_scale = create_colour_scale(3, max_n_species);

//     vector<double> colour_0{ 68, 13, 86 };
//     vector<double> colour_1{ 32, 143, 141 };
//     vector<double> colour_2{ 247, 229, 31 };
//     vector<vector<double>> colour(3, vector<double>(max_n_species, 0));
    
    
//     for(int i_number = 0; i_number < max_n_species; i_number++){
//         colour[0][i_number] = colour_scale[0][i_number]*colour_0[0]
//                             + colour_scale[1][i_number]*colour_1[0]
//                             + colour_scale[2][i_number]*colour_2[0];

//         colour[1][i_number] = colour_scale[0][i_number]*colour_0[1]
//                             + colour_scale[1][i_number]*colour_1[1]
//                             + colour_scale[2][i_number]*colour_2[1];

//         colour[2][i_number] = colour_scale[0][i_number]*colour_0[2]
//                             + colour_scale[1][i_number]*colour_1[2]
//                             + colour_scale[2][i_number]*colour_2[2];
//     }
    



//     colorscale_R = 3*255.0 / log(max_n_species+1);
//     colorscale_G = 255.0 / log(max_n_species+1);
//     colorscale_B = 2*255.0 / log(max_n_species+1);


    
//     int n_message = 4 + model_ptr->n_reactions + model_ptr->n_species + 2;
//     string state_value, messages[4 + MAX_REACTIONS + MAX_SPECIES + 1];
    

//     sf::RectangleShape rectangle(sf::Vector2f(s_size, s_size));
//     sf::RectangleShape legend_box(sf::Vector2f(l_size, l_size));

//     sf::Text text;
//     text.setCharacterSize(text_size);
//     text.setFillColor(sf::Color::White);

//     sf::Font font;
//     if (!font.loadFromFile("Fonts/arial.ttf")) {
//         cout << "Could not load font..." << endl;
//     }

//     text.setFont(font);
// 	text.setOutlineColor(sf::Color::Black);	        
// 	text.setOutlineThickness(2);	

// 	sf::Vector2i localPosition; 
	
//     while (window.isOpen() ){

//         window.clear();

//         // Drawing text

//         messages[0] = "Elapsed clock time:" + to_string(clock.getElapsedTime().asSeconds());   // Elapsed time
//         messages[1] = "Elapsed simulation time:" + to_string(model_ptr->t);                    // Elapsed simulation time
//         messages[2] = "Iteration number:" + to_string(model_ptr->i_iter);                      // Elapsed iterations
//         messages[3] = "";

//         // Counts of all types of reactions
//         for(int i_reaction = 0; i_reaction < model_ptr->n_reactions; i_reaction++){
//             reaction_ptr = &(model_ptr->reactions[i_reaction]);
//             messages[4+i_reaction] = reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count);
//         }
//         //todo find how to change value per cell and not the same over and over
//         // count of total amount of nanostars bound
//         for(int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {
            
//             cell_ptr = &(model_ptr->grid[i_cell]);
//             int bound_nanostars = 0;
//             // // calculate total number of bound nanostars
//             // bound_nanostars = 25000000 - model_ptr->n_species-5
//             //for every bound species add amount to total
//             for(int i_species = 0; i_species < model_ptr->n_species; i_species++){
//                 //if nanostar bound add to total of bound nanostars
//                 if(model_ptr->species[i_species].find('1') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('2') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('3') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('4') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('5') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('6') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('7') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('8') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 } else if(model_ptr->species[i_species].find('9') != std::string::npos){
//                     bound_nanostars += cell_ptr->state[i_species];
//                 }
//             }

//             //calculate coverage
//             int max_coverage = 10000; //just put a number as when you take the logarithm this constant wont matter ( power )
//             double coverage = static_cast<double>(bound_nanostars)/static_cast<double>(max_coverage);

//             //store nanostars_tot and coverage (nanostars_tot is the third to last species)
//             cell_ptr->state[model_ptr->n_species-2] = bound_nanostars;
//             cell_ptr->state[model_ptr->n_species-1] = coverage;

//             // for(int i_species = 0; i_species < model_ptr->n_species; i_species++){
//             //     if(model_ptr->species[i_species].find('nanostars') != std::string::npos){
//             //         // cell_ptr->state[i_species] = bound_nanostars;
//             //         cout<<to_string(cell_ptr->state[i_species]) +": "+to_string(bound_nanostars)+"\n";
//             //     }
//             // }
//         }

//         // //alpha causes too many problems
//         // for(int i_cell = 0; i_cell < model_ptr->n_cells-1; i_cell++) {

//         //     cell_ptr = &(model_ptr->grid[i_cell]);
//         //     double theta1 = cell_ptr->state[model_ptr->n_species-2];
//         //     double dN1 = cell_ptr->state[model_ptr->n_species-4];

//         //     cell_ptr = &(model_ptr->grid[i_cell+1]);
//         //     double theta2 = cell_ptr->state[model_ptr->n_species-2];
//         //     double dN2 = cell_ptr->state[model_ptr->n_species-4];

//         //     //calculate dtheta
//         //     double dtheta = log(theta2) - log(theta1);

//         //     //calculate dNr
//         //     double dN = log(static_cast<double>(dN2)) - log(static_cast<double>(dN1));

//         //     //calculate and store alpha
//         //     cell_ptr->state[model_ptr->n_species-1] = dtheta/dN;
//         // }

//         // // Counts of all types of diffusion
//         // for(int i_species = 0; i_species < model_ptr->n_species; i_species++){
//         //     reaction_ptr = &(model_ptr->diffusion_reactions[i_species]);
//         //     messages[4+model_ptr->n_reactions+i_species] = reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count);
//         // }

//         //     reaction_ptr = &(model_ptr->diffusion_reactions[i_species]);
//         //     messages[4+model_ptr->n_reactions+i_species] = reaction_ptr->reaction_name + " count: " + to_string(reaction_ptr->count);
//         // }

//         messages[4+model_ptr->n_reactions] = "";
//         messages[5+model_ptr->n_reactions] = "Showing species " + model_ptr->species[plot_species%model_ptr->n_species];

//         // Drawing messages
//         for(int i_message = 0; i_message < n_message; i_message++){
//             text.setString(messages[i_message]);
//             text.setPosition(border + height + 0.5*border + (l_size + 2*s_border) + 0.5*border, border + i_message * line_size);
//             window.draw(text);
//         }

//         // TODO: read new code
//         //create color legend
//         x0 = border + height + 0.5*border; 

// 		for (int i_legend = 0; i_legend < N_legend; i_legend++) {
		
// 			y0 = border + s_border + i_legend * (l_size + 2*s_border); 
		
// 	    	legend_box.setPosition(x0,y0);
// 	    	text.setPosition(x0,y0);
	    	
// 	    	factor = pow(max_n_species,i_legend / double(N_legend-1));
// 	    	text.setString(to_string(int(factor)));
// 			// legend_box.setFillColor(sf::Color(int(colorscale_R*log(factor+1))%256, int(colorscale_G * log(factor+1))%256, int(colorscale_B * log(factor+1))%256));
//             // legend_box.setFillColor(sf::Color(min(int(colorscale_R*log(factor+1)), 255), min(int(colorscale_G * log(factor+1)),255), min(int(colorscale_B * log(factor+1)), 255)));
//             legend_box.setFillColor(sf::Color(colour[0][int(factor)], colour[1][int(factor)], colour[2][int(factor)]));
// 			window.draw(legend_box);
// 			window.draw(text);
// 		}

//         // Getting position of mouse cursor.
// 		localPosition = sf::Mouse::getPosition(window);
//         // END TODO: read new code

//         // Drawing cells
//         for(int i_cell = 0; i_cell < model_ptr->n_cells; i_cell++) {

//             cell_ptr = &(model_ptr->grid[i_cell]);
            
//             // Top right coordinate of the current cell.
// 			x0 = border + s_border + scale*(cell_ptr->x - model_ptr->min_x);
// 			y0 = border + s_border + scale*(cell_ptr->y - model_ptr->min_y);
            
//             rectangle.setPosition(x0,y0);
            
//             // cout << "Position top: " << border + s_border + scale*(cell_ptr->y - model_ptr->min_y) << endl;
//             // cout << "Position bottom: " << border + s_border+ scale*(cell_ptr->y - model_ptr->min_y) + s_size<< endl;
//             // cin.get();

//             // TODO: read new code
//             count_species = cell_ptr->state[plot_species%model_ptr->n_species];
//             rectangle.setFillColor(sf::Color(colour[0][int(count_species)], colour[1][int(count_species)], colour[2][int(count_species)]));
//             // rectangle.setFillColor(sf::Color(int(colorscale_R * count_species)%256,int(colorscale_G * count_species)%256,int(colorscale_B * count_species)%256));
//             // END TODO: read new code

//             // state_value = to_string(c);
//             // text.setString(state_value);
//             // text.setPosition(scale*50.f + scale*35.f * cell_ptr->x, scale*50.f + scale*35.f * cell_ptr->y);
            
//             window.draw(rectangle);

//             // TODO: read new code
//             // If mouse in the current cell, display current contents of that cell.
//             if (localPosition.x > x0 and localPosition.x < (x0 + s_size) and localPosition.y > y0 and localPosition.y < (y0 + s_size)) {
//                 LocalComposition = "";
//                 //print #molecules as an integer
//                 // for(int i_species = 0; i_species < model_ptr->n_species-2; i_species++){
// 	            //     LocalComposition += model_ptr->species[i_species] + " " + to_string(static_cast<int>(cell_ptr->state[i_species])) + "\n";
// 	            // }
//                 LocalComposition += model_ptr->species[model_ptr->n_species-2] + " " + to_string(cell_ptr->state[model_ptr->n_species-2]) + "\n";
//                 LocalComposition += model_ptr->species[model_ptr->n_species-1] + " " + to_string(cell_ptr->state[model_ptr->n_species-1]) + "\n";
//                 LocalComposition += "time: " + to_string(model_ptr->t)+ "\n";

// 				text.setString(LocalComposition);
// 				text.setPosition(localPosition.x+30,localPosition.y-30);        
// 	        }
//             // END TODO: read new code
//         }

//         window.draw(text);
//         // TODO: remove time_h clock, use internal sf clock.
//         // For switching which species is currently rendered.
//         // It checks whether the left mouse button has been clicked. If it is clicked, and time has advanced more than 200 ms,
//         // then make change to the plotted species and reset click time. Ensures held down click does not advance species.
// 		if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
// 			Current = clock.getElapsedTime(); 
// 			if (Current.asMilliseconds() - LeftClick.asMilliseconds() > 200) {
// 				plot_species ++;//should be ++
// 				LeftClick = clock.getElapsedTime(); 
// 			}
// 		}
// 		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
// 			Current = clock.getElapsedTime(); 
// 			if (Current.asMilliseconds() - RightClick.asMilliseconds() > 200) {
// 				plot_species --;//should be --
// 				RightClick = clock.getElapsedTime(); 
// 			}
// 		}
    	

//         sf::Event event;
//         while (window.pollEvent(event))
//         {
//             if (event.type == sf::Event::Closed)
//                 window.close();
//         }

//         window.display();
//     }
// }


void make_model(Model* model_ptr){

    time_t t_start, t_end;
    int clock_time[3];

    model_ptr->create_cartesian_grid(1, 13, 1, 1, 13, false);//last two ints determine size of the grid //Bookmark

    model_ptr->load_model("6 arms konsol 1.00e-15 konsur 5.00e-04_model.txt"); //! This is where the model.txt file is specified
    model_ptr->load_initial_conditions("typeAB_setup.txt");

    time(&t_start);

    //TODO: make this a function and make sure scale can easily be changed 
    //set receptors in cells
    float low_power = 2;
    float high_power = 4.5;
    float n_powers = high_power-low_power;
    float n_recept_max = static_cast<int>(pow(10, high_power));
    long long int n_recept_array [13] = {};
	long long int n_b_array[13] ={};
    //some pointer n_cells
    for(int ii = 0; ii < model_ptr->n_cells; ii++){
        double receptorpower = ((static_cast<double>(ii*n_powers))/12) + low_power;
        n_recept_array[ii] = static_cast<long long int>(pow(10, receptorpower));
		// n_b_array[ii] = static_cast<long long int>(pow(10, receptorpower));
		n_b_array[ii] = static_cast<int>(0);        
        //n_b_array[25-1-ii] = static_cast<int>(pow(10, receptorpower));
    }
    
    //Boqi made
    cout << "number of receptors in each cell is: ";

    int n = sizeof(n_recept_array) / sizeof(n_recept_array[0]);
    for (int i = 0; i < n; i++)
        cout << n_recept_array[i] << endl;

    //end of what Boqi made

    for(int ii = 0; ii < model_ptr->n_cells; ii++){
        model_ptr->grid[ii].set_species_number(0, n_recept_array[ii]); // here species number 0 points to the first species defined in model.txt (recept1) is set to concentration n+receptors
        model_ptr->grid[ii].set_species_number(model_ptr->n_species-4, n_recept_array[ii]); // Here the total number of recept1 is set to n_receptors 
        model_ptr->grid[ii].set_species_number(1, n_b_array[ii]);// here species number 0 points to the first species defined in model.txt (recept1) is set to concentration n+receptors
        model_ptr->grid[ii].set_species_number(model_ptr->n_species-3, n_b_array[ii]);
    }

        
    //BOQI made
    // Setting file name to include project and cell index.
    string folder_name;
    folder_name = "./Data/";
    string file_name = folder_name + string("parameters") + ".dat";

	ofstream datafile;
    datafile.open(file_name);
    datafile << "low power" << setw(20) << low_power << endl;
    datafile << "high power" << setw(20) << high_power << endl;
    datafile.close();

    // 100 million iterations
    //model_ptr->run(100000000, 100000, false);
    //put for loop here for kc/koff
    //put for loop here to do parameter sweep over kc/ko
    // long run
    //model_ptr->run(1000000000, 100000, false);
    // model_ptr->run(50000000000, 10000, false, folder_name);
    model_ptr->run(50000000000, 10000, false, folder_name);
    
    // medium long run
    //model_ptr->run(100000000, 100000, false);
    // short run
    // model_ptr->run(10000000, 10000, false);


    time(&t_end);

    // model_ptr->grid[13].set_species_number(0, 10000);
    // model_ptr->grid[0].set_species_number(1, 130);
    // model_ptr->grid[1].set_species_number(0, 500);
    


    // TODO: also check events in event creation loop.
    // Checking reactions

    // Reaction* reaction_ptr;
    // for(int i = 0; i < model_ptr->n_reactions; i++){
    //     reaction_ptr = &(model_ptr->reactions[i]);

    //     cout << reaction_ptr->reaction_name << endl;

    //     for(int j = 0; j < reaction_ptr->n_reaction_species;j++){
    //         cout << "Species "<< reaction_ptr->species[j] << " changes by " <<
    //         reaction_ptr->state_change[j] << " with rate " << reaction_ptr->k << endl;
    //     }
    // }



    // Saving hours, minutes, seconds.
    clock_time[0] = int(difftime(t_end,t_start)/3600.0);
    clock_time[1] = int(difftime(t_end,t_start))%3600/60;
    clock_time[2] = int(difftime(t_end,t_start))%60;

    
    // Checking Reactions
    // for(int j = 0; j<3; j++){
    // for(int i = 0; i<model.reactions[j].n_reaction_species; i++){
    //     cout << "species "<< model.reactions[j].species[i] << " changes by " << model.reactions[j].state_change[i] << "\n";
    // }
    // }

    // Checking neighbours
    // for(int i = 0; i<model_ptr->n_cells; i++){

    //     cout << "cell " << i << " neighbours \n";

    //     for(int j = 0; j<model_ptr->grid[i].n_neighbours; j++){
            
    //         cout << "Neighbour: " << model_ptr->grid[i].neighbours[j] << endl;
    //         cout << "Distance: " << model_ptr->grid[i].neighbour_dist[j] << endl;

    //     }
    //     cout << endl;
    // }
}

int main() {

    // Put model in heap memory, so that it can be larger. 
    Model* model_ptr = new Model;
    //sf::Thread thread(&make_model, model_ptr); //Mac does not allow the opening of windows that are not in the main thread, for ... reasons, so the simulations occurs in a separate thread
    //thread.launch();
    //render(model_ptr);

    make_model(model_ptr);
	
    // Model* model_ptr = new Model;
    // make_model(model_ptr);;
    
    cout << "Press any key to continue";
    cin.get();
    return 0;
}
