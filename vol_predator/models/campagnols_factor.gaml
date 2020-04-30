/**
 *  Voles model
 *  Author: nicolas
 *  Description: 
 *  On travail uniquement sur les femelles, pour avoir des résultats globaux, il faut multiplier par 2. Les prédateurs mange donc, en femelle, la moitier de leurs besoins alimentaires. 
 * 
 * 
 */

model campagnols

global {
	// Space projection (boundaries of the space)
	file emprise <- file("../includes/emprise1.shp");
	// Elevation model of the space
	file space_grid_file <- file("../includes/DEM_Romanche_LIIres50.asc");
	//  studied space
	file plot_shape <- file("../includes/Mailles_2014_reduit.shp");
	// Source place of vol colonisation
	file starting_shape <- file("../includes/historique_1998.shp");
	// predefined landuse color
	map colors <- map([1::rgb([178, 180, 176]), 2::rgb([246, 111, 0]), 3::rgb([107, 0, 0]), 4::rgb([249, 0, 255]), 5::rgb([144, 96, 22]), 6::rgb([255, 255, 86]), 7::rgb([19, 114, 38]), 8::rgb("black"), 9::rgb([107, 94, 255]), 10::rgb([43, 255, 255])]);
	
	geometry shape <- envelope(emprise);
	list<cells> simulated_plots <- [];


///////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//poids des campagnols
	float young_weight <- 40#gram;
	float juvenile_weight <- 80#gram;
	float adult_weight <- 80#gram;
	// entre 2 et 11 kg par KM2/J
	
	float qt_female_to_eat <- 500#kilo/(#km*#km)/#day;    ////*30*/ (5 * (80#gram))/(#km*#km);
	//float vole_eat_qte <-30;
	//float vole_eat_surfacique_mass <- vole_eat_qte*80#gram/(100#m*100#m);
	float adult_eat_proportion <- 1.0;//0.5;
	float juvenile_eat_proportion <- 0.0;//0.5;
	float young_eat_proportion <- 1 - adult_eat_proportion - juvenile_eat_proportion;

	float buffer_size <- 200#m;





///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int simnum<-0;
	
	// maximum age of a younger vol
	int younger_compartment_max_age <- int(18#day) ;
	// maximum age of a juvenile vol
	int juvenile_compartment_max_age <- int((8*7)#day) ;
	//grow rate of younger vol from adult and juvenile vol. (individual/date/individual)
	float adult_grow_rate <- 0.0165/#day; 
	
	/**
	 * Death rates.
	 */
	//Death rate of younger vol.  (individual/date/individual)
	float younger_death_rate <- 0.0;
	//Death rate of juvenile vol.  (individual/date/individual)
	float juvenile_death_rate <-0.0;
	//Death rate of adult vol.  (individual/date/individual)
	float adult_death_rate <- 0.0; 
	//vol density saturation threshold
	float dispersion_threshold <- 0.2;
	
	//Structural capacity of the space in indiividuals / Ha
	float default_vole_carrying_capacity <- 500/(100#m*100#m);

	//rate for selecting vol behavior (short range motion vs long range motion)
	float crazy_vole_rate <- 0.1;
	//Long range moving distance : minima and maxima
	float min_distance_crazy_vole <- 3#km;
	float max_distance_crazy_vole <- 6#km;
	// moving  behaviour probability (probability of moving to an equal or lower elevation)  	
	float coeff_go_down <- 2/3;
	
	//reproducing season (from the begining of april to the end of october)
	int open_season <- 4;
	int close_season <- 10;
	
	/**
	 * starting date of the simulation
	 */
	date dstart <- date([1998, 4,1,0,0,0]);
	//date dstart <- date([1998, 1,1,0,0,0]);
	
	/**
	 * internal variables of the model 
	 */
	float born_accumulator <- 0.0;
	float death_accumulator <-0.0;
	
	/**
	 * variables containing outputs of the model
	 */
	float pop_total <- 1.0 update: (sum(cells collect(each.total_population)));
	float pop_adult <- 0.0 update: sum(cells collect(each.adult_accumulator));
	float pop_juvenile <- 0.0 update: sum(cells collect(each.juvenile_accumulator));
	float pop_younger <- 0.0 update:sum(cells collect(each.younger_accumulator));
	float colonisation_speed <- 0.0;
	float discovered_area <- 0.0;
	float discovered_area_year <- 0.0;
	float densite_monitor <- 0.0 update:mean(cells collect(((each.total_population/each.shape.area)/default_vole_carrying_capacity)*100));
	float densite_occuped_monitor <- 0.0 update:mean(cells where(each.is_occupated) collect(((each.total_population/each.shape.area)/default_vole_carrying_capacity)*100));
	float population_young_monitor <- 0.0 update:sum(cells collect(each.younger_accumulator));
	float population_juvenile_monitor <- 0.0 update:sum(cells collect(each.juvenile_accumulator));
	float population_adult_monitor <- 0.0 update:sum(cells collect(each.adult_accumulator));
	float max_cell_population <- 0.0 update: max(cells collect(each.total_population));
	int long_distance_dispersing_voles <- 0 update:(voles count(each.crazy_voles));
	int short_distance_dispersing_voles <- 0 update:(voles count(!each.crazy_voles));
	int dispersing_voles -> {short_distance_dispersing_voles + long_distance_dispersing_voles };
	
	
	init
	{
		// one simulation step correspond to one day in the real life
		step <- 1#day;
		starting_date <- dstart;
		/**
		 * Creation of the space
		 */
		create plots from:plot_shape;
		loop pp over:plots
		{
			 ask cells overlapping(pp)
			{
				carring_capacity <- self.shape.area*default_vole_carrying_capacity;
				color <- #green;
			}
		}
		simulated_plots <- cells where(each.carring_capacity != 0);
		
		/**
		 * creation of vole population outbreak starting point
		 */
		create voles_start_point from: starting_shape 
		{
			adult_init <- 100;
			younger_init <- 0;
			juvenile_init <- 0;
			
			write "adulte " + adult_init;
			write "juvenile " + juvenile_init;
			write "younger" + younger_init;
			write "grow " + adult_grow_rate;
		}
		
		ask voles_start_point
		{			
			ask cells overlapping(self)
			{
				color <- #red; 
				adult_accumulator <- float(myself.adult_init);
				int i <- 0;
				young_born_day_acc<- [];				
			}
		}

		/**
		 *  initialise the elevation model
		 */
 		ask mnt
		{
				associated_cell <- one_of(cells overlapping self);
				associated_cell.my_mnt_cell <- self;
				neighbours_mnt_cell <- (self neighbors_at 1 ) ;
				int down <- neighbours_mnt_cell count(each.height <=self.height);
				sum_weight <- coeff_go_down * down + (1-coeff_go_down) * (length(neighbours_mnt_cell) - down);
		}
		
		ask simulated_plots
		{
			my_far_cell <- simulated_plots overlapping (max_distance_crazy_vole around circle(min_distance_crazy_vole));
		}
		
		
	
		create death_factor number:1
		{
			location <- one_of(simulated_plots).location;
			my_color <- rgb(rnd(255),rnd(255),rnd(255));
		}
		
		ask simulated_plots
		{
			cells c_cell <- self;
			ask death_factor with_min_of(each.location distance_to c_cell.location  ) 
			{
				do add_to_knowledge(c_cell);
			}	
		}
		
		ask death_factor{
			do update_boundaries;
			qt_vole_to_eat <-1#gram ; //500#kilo; //(length(knowledge.keys)*first(knowledge.keys).shape.area)* qt_female_to_eat; //100#gram;//
			write "superficie qte "+ qt_female_to_eat;
			write "superficie surf  "+ (length(knowledge.keys)*first(knowledge.keys).shape.area);
			write "superficie "+ ((length(knowledge.keys)*first(knowledge.keys).shape.area))+ " " + qt_vole_to_eat;
			self.adult_eat_proportion <- world.adult_eat_proportion;
			self.juvenile_eat_proportion <- world.juvenile_eat_proportion;
			self.young_eat_proportion <- world.young_eat_proportion;
	
		}
 
		
	}
	
	reflex update_bounds
	{
		ask death_factor{
			do update_boundaries;
		}
	}
	
	/**
	 * stop simulation when all vols are destroyed
	 */
	reflex no_vol when: false and population_young_monitor + population_juvenile_monitor + population_adult_monitor <= 0
	{
		do halt;
	}
	
	/**
	 * action to do at the begining of a year
	 */
	reflex new_year when:  dstart.month = current_date.month and dstart.day = current_date.day and current_date.year != dstart.year
	{
		float area_cell <- one_of(cells).shape.area;
		float discovered_cell <- cells count (each.is_occupated) * area_cell;
		discovered_area_year <- discovered_cell - discovered_area;
		discovered_area <- discovered_cell;
	} 
	
	/**
	 * action to do at the begining of a semester
	 */
	reflex newsemester when: current_date.day=1 and ( current_date.month = 6 or current_date.month = 12)
	{
	save simulated_plots with:[younger_accumulator::"younger",juvenile_accumulator::"juvenile",adult_accumulator::"adult"] to:"outs/cells_"+simnum+"_"+cycle+".shp" type:"shp";
	if(length(voles)>0) {save voles to:"outs/voles_"+simnum+"_"+cycle+".shp" type:"shp";}
	} 
	
}

/**
 * specie to define vol pululation starting point
 */
species voles_start_point
{
	int younger_init <- 0;
	int juvenile_init <- 0;
	int adult_init <- 0;
	aspect default
	{
		draw circle(30#m) color:#blue;
	}
}

species zone_information schedules:[]
{
	
	date last_appliance; //last date when it eat in the cell
	date date_to_forget; //date when the knowledge is forgoten
	float number_of_young; //the number of animals it has seen
	float number_of_juvenile;
	float number_of_adult;
	float number_of_voles -> {number_of_young + number_of_juvenile + number_of_adult };
	int last_update_cycle <- -1;
	cells associated_cell;
	
	bool is_boundary <- true;
	
	action update
	{
		last_appliance <- current_date;
		date_to_forget <- last_appliance + gauss(30#day,5#day);
		number_of_young <- associated_cell.younger_accumulator;
		number_of_juvenile <- associated_cell.juvenile_accumulator;
		number_of_adult <- associated_cell.adult_accumulator;
		last_update_cycle <- cycle;
	}
	
	float evaluate 
	{
		return 0;
	}
	
	aspect base
	{
		if(cycle-1 = last_update_cycle)
		{
			draw associated_cell.shape color:#blue;
		}
	}
	 
}
// model predators of the system. this species should be inherited to derive specific predators (generalist, specialist ...)
species death_factor
{
	float qt_vole_to_eat;  //en gramme par jour
	
	float adult_eat_proportion <- 0.0;
	float juvenile_eat_proportion <- 0.0;
	float young_eat_proportion <-0.0;
	
	list<zone_information> used_knowledge <-[];
	
	//proportion de chaque quan
	
	int travel_to_eat_size <- 10; //0000;
//	float population;
	map<cells,zone_information> knowledge <- [];	
	rgb my_color ;
	
	geometry my_research_area <- location;
	
	reflex reduce_zone when:true
	{
		list<zone_information> to_forget <- knowledge.values where(each.date_to_forget<=current_date and each.is_boundary = true  );
		
		// rechercher les nouvelles cellules boundaries dans les boundaries...
		
		list<zone_information> to_forget_ordered <- knowledge.values sort_by (each.date_to_forget);
		write "to remove "+ length(to_forget)+" "+ current_date +" "+first(to_forget_ordered).date_to_forget;	
		loop info over:to_forget{
			if(info.is_boundary = true)
			{
				list<cells> ml <-knowledge.keys inter (info.associated_cell.neighbors);
				loop cell over: ml
				{
					((knowledge at cell).is_boundary) <- true;
				}
			}
			remove key:info.associated_cell from:knowledge;
		}
	}
	
	reflex discover when: (!empty(used_knowledge) )
	{
		write "discover "+length(used_knowledge)+ " "+length(knowledge);
		//il découvre des cellules (nouvelles ou de sa base de connaissances) quand il n'a plus assez à manger
		// trace le polygone P convexe entre les n cellules où il s'est nourri.
		// on fait une homothétie H sur ce polygone P avec un rapport calculé pour que la surface de ce nouveau polygone soit égal à Alpha * la surface du polygone de voronoi de départ. Alpha paramètre du modèle. 
		// dans le polygone H on actualise toute les connaissances et construit de nouvelles connaissances. 
		//	used_knowledge <- knowledge;
		list<point> all_location <- used_knowledge collect(each.associated_cell.location);
		
		my_research_area <- convex_hull(polygon(all_location));
		// il met a jour tout ces nouvelles cellules
		my_research_area <-  union( [my_research_area, ( buffer_size around my_research_area)]);
		list<cells> all_cell <- simulated_plots overlapping my_research_area;
		loop a over:all_cell
		{
			do add_to_knowledge(a);
		}
	
		//j'agrandis mon territoire si et seulement si je suis en disette cad avoir moins de campagnols sur ma base de connaissance par rapport à ce dont j'ai besoin
		used_knowledge <- [];
	}
	
	reflex eat when:cycle>1
	{
		// OK -- les prédateurs sont qualifiés par une qte de campagnol à manger chaque jour.
		// OK -- ils vont classer les cellules en fonction du nombre de campagnols disponibles
		list<zone_information> sorted_cells <- knowledge.values sort(each.number_of_voles);
		int nb_cell_moved <-0;
		int nb_cell_moved_start <-length(sorted_cells);
		float t_to_eat <- qt_vole_to_eat;
		
		//loop while:t_to_eat > 0 and nb_cell_moved < nb_cell_moved_start
		//{
			int knowledge_size <- length(sorted_cells);
			list<zone_information> first_elements <- copy_between(sorted_cells,knowledge_size-travel_to_eat_size,knowledge_size);
			// ils mangent dans les n premières cellules au prorata du nombre de campagnols dans chaque cellule
			
			write "taille elements " + knowledge_size+" " +length(first_elements);
			write "data "+ first_elements collect(each.associated_cell.total_mass);
			float available_meat <- sum(first_elements collect(each.associated_cell.total_mass));
			list<float> killed <- [0.0,0.0,0.0];
			list<float> local_killed <- [0.0,0.0,0.0];
			nb_cell_moved <- nb_cell_moved +length(first_elements);
			if(available_meat > 0){
				ask first_elements
				{
					// les bébés sont ils mangé ? oui non, un paramètre à ajouter...
					float to_eat  <- t_to_eat * (self.associated_cell.total_mass/available_meat);
					ask self.associated_cell {
						//correction --- parler en individus
					//	if(available_meat > to_eat)
					//	{
							local_killed <- kill(to_eat*adult_eat_proportion,to_eat*juvenile_eat_proportion, to_eat*young_eat_proportion);
					//	}	
					//	else
					//	{
					//		local_killed <-  kill_all(adult_eat_proportion!=0.0, juvenile_eat_proportion!=0.0,young_eat_proportion!=0.0);
					//	}
						//t_to_eat <- t_to_eat - sum(local_killed);
					}
					// S'il n'a plus assez a manger dans les n premières cellules (1ere iteration), il doit actualiser sa base de connaissance par une exploration
					// il va chercher dans d'autres cellules de sa base de connaissances ou actualiser sa base de connaissances.
					
				}
			}
			remove all:first_elements from: sorted_cells;
			add all:first_elements to:used_knowledge;
		//}
		
//		write "nb cell " + nb_cell_moved +  " "+ nb_cell_moved_start;
//		write  self.name +" want to eat "+t_to_eat+ "("+qt_vole_to_eat+")";
			
	}
	
	
	
	
	bool is_boundary_cell(cells cell_to_test)
	{
		list<cells> tmp <- cell_to_test.neighbors;
		bool res <- knowledge.keys contains_all(tmp);
		return !res;
	}
	
	action add_to_knowledge(cells my_cell)
	{
		int aa <- length(knowledge.keys);
		zone_information my_know <- nil;
		if(knowledge.keys contains(my_cell)) {
			zone_information my_know <- knowledge[my_cell];
			ask my_know
			{
				do update;
			}
		}
		else {
			write "add zone ";
			create zone_information number:1
			{
					associated_cell <- my_cell;
					self.shape <- my_cell.shape;
					self.location <- my_cell.location;
					
					do update;
					my_know <- self;
			}
			//add my_know at:my_cell to:knowledge;
			
			knowledge <-knowledge + [my_cell::my_know];
			int bb <- length(knowledge.keys);
		}
	}
	
	action update_boundaries {
		loop cell over:knowledge.keys
		{
			(knowledge at cell).is_boundary <- is_boundary_cell(cell);
		}
	}
	
	action add_zone(cells mc)
	{
		if(knowledge[mc] = nil)
		{
			create zone_information with:[number_of_young::0,number_of_juvenile::0,number_of_adult::0, last_appliance::current_date, associated_cell::mc] returns:agt;
			add first(agt) at:mc to:knowledge;
		}
	}
	
	
	aspect base
	{
		cells mc;
		loop mc over:knowledge.keys
		{
			draw mc.shape color:(knowledge at mc).is_boundary?#black:my_color ;	
		}
		
		draw circle(100#m) color:#red;
	}
	
	aspect discovering_aspect
	{
		draw my_research_area color:my_color;
	}
	
}


/**
 * Vole species, it contains agent characteristics and agent behaviour (move, die)
 */
species voles 
{
	//age of the juvenile agent
	int juvenile_age_id <- 0;
	//does the voles follow long distance behaviour or short distance distance behaviour
	bool crazy_voles<-true;
	//
	mnt current_place <- nil;
	cells my_current_cell;
	
	//long distance moving behaviour
	reflex gogo_crazy when: crazy_voles
	{
		my_current_cell <- one_of(my_current_cell.my_far_cell); //one_of(simulated_plots overlapping (6#km around circle(3#km)));
		location <- my_current_cell.location;
	}
	
	//short distance moving behaviour
	reflex gogo when: !crazy_voles
	{
		//list<mnt> places <- current_place.neighbours_mnt_cell; //neighbours_at 1 where(each.height <=current_place.height);
		mnt plc <- choose_cell(current_place,rnd(1000)/1000,coeff_go_down,current_place.sum_weight);    //length(places)>0?one_of(places):current_place;
		current_place <- plc;
		location <- current_place.location;
		my_current_cell <- current_place.associated_cell;// one_of(cells overlapping current_place);
	}
	
	//determine the next cell to move according to the current spatial situation
	mnt choose_cell(mnt mplace, float rndVal, float coeff_down, float maxVal)
	{
		list<mnt> places <- mplace.neighbours_mnt_cell;//(self neighbours_at 1 ); //where(each.height <=self.height )
		float val;
		mnt res;
		loop lc_place over:places
		{
			val <- val + (mplace.height>=lc_place.height ? coeff_down:(1-coeff_down))/ maxVal;
			if(val>=rndVal)
			{
				res <- lc_place;
				break;
			}	
		}
		return res;
	}
	
	// vulnaribility to die on the way to the next cell (predation, ...)
	reflex killed when:flip(0.5)
	{
		do die;
	}
	
	// The vole find a place to settle down, it is integrated into age category
	reflex individual_to_population when: my_current_cell!=nil and my_current_cell.carring_capacity!= 0 and my_current_cell.total_population/my_current_cell.carring_capacity<dispersion_threshold
	{
		put my_current_cell.juvenile_day_acc at juvenile_age_id + 1 at: juvenile_age_id in: my_current_cell.juvenile_day_acc;  
		do die;
	}
	
	// Just for display 
	aspect default
	{
		draw circle(20#m) color:#blue;
	}
}

/**
 * DEM grid 
 */
grid mnt   file: space_grid_file schedules:[] {
	float height <- grid_value;
	cells associated_cell <- nil;
	list< mnt> neighbours_mnt_cell  <- nil;
	float sum_weight <- 0.0;	
}

/**
 * population grid, each cell is an agent (cell agent) that control local population
 */
grid cells cell_width:100#m cell_height:100#m schedules:  simulated_plots
{
	// local carring capacity
	float carring_capacity <- 0.0;

	list<cells> my_far_cell<-[];
	mnt my_mnt_cell;
	//daily accumulator by actual age
	list<float> young_born_day_acc <-[];
	list<float> juvenile_day_acc <-[];

	//daily accumulator by age class
	float younger_accumulator -> {sum(young_born_day_acc)} min:0.0;
	float juvenile_accumulator -> {sum(juvenile_day_acc)} min:0.0;
	float adult_accumulator <- 0.0 min:0.0;

	//cell monitors
	float total_population -> {sum(young_born_day_acc) + sum(juvenile_day_acc) +adult_accumulator };
	float total_mass -> {sum(young_born_day_acc)*young_weight + sum(juvenile_day_acc)*juvenile_weight  +adult_accumulator*adult_weight  };
	
	float adult_rate -> {adult_accumulator = 0? 0:adult_accumulator / total_population};

	rgb color <- #white update: carring_capacity=0?#white:(total_population!=0?#red:#green);
	bool is_occupated -> {total_population>0};
	
	
	aspect predating_zone
	{
			draw shape color:#white;
	}
	
	/**
	 * for the display
	 */
	aspect grille
	{
		if(carring_capacity != 0)
		{
			draw shape color:#white;			
		}
	}
	
	/**
	 * FIFO accumulator managment - get the older out
	 */
	float pull_data(list<float> mList,int da)
	{
		float data <- 0.0;
		if(length(mList)>=(da/step))
		{
			data <- mList[0];
			remove index:0 from:mList;
		}
		return data;
	}
	
	list<float> kill_all(bool adult, bool juvenile, bool young)
	{
		float res_qte_younger <- 0.0;
		float res_qte_adult <- 0.0;
		float res_qte_juve <- 0.0;
		if(juvenile)
		{
			int i <-0;
			loop while:i<length(juvenile_day_acc)
			{
				res_qte_juve <- res_qte_juve + juvenile_day_acc[i];
				juvenile_day_acc[i] <- 0;
				i <- i +1;
			}	
		}
		if(young)
		{
			int i <-0;
			loop while:i<length(young_born_day_acc)
			{
				res_qte_younger <- res_qte_younger + young_born_day_acc[i];
				young_born_day_acc[i] <- 0;
				i <- i +1;
			}	
		}
		if(adult)
		{
			res_qte_adult <- adult_accumulator;
			adult_accumulator <- 0.0;	
		}
		list<float> to_ret<- convert_ind2weight(res_qte_adult, res_qte_juve, res_qte_younger);
		return to_ret;
	}
	
	// convertit des masses en individus
	list<float> convert_weight2ind(float adult,float juve, float young) {
		float i_adult <- adult/adult_weight;
		float i_juvenile <- juve / juvenile_weight;
		float i_young <- young /young_weight;
		return list([i_adult,i_juvenile,i_young]);
	}
	// convertit des individus en masse
	list<float> convert_ind2weight(float adult,float juve, float young) {
		float i_adult <- adult*adult_weight;
		float i_juvenile <- juve * juvenile_weight;
		float i_young <- young *young_weight;
		return list([i_adult,i_juvenile,i_young]);
	}
	
	// paramètre en gramme et retourne une masse en gramme...
	list<float> kill(float g_adult,float g_juve, float g_young) {
		list<float> ml <- convert_weight2ind( g_adult, g_juve,  g_young);
		
		
		
		float adult <- min([ml[0],adult_accumulator]);
		float juve <- min([ml[1],juvenile_accumulator]);
		float young <- min([ml[2],younger_accumulator]);
		
		float res_qte_adult <- adult;
		float res_qte_juve <- 0.0;
		float res_qte_young <- 0.0;

		adult_accumulator <- adult_accumulator - adult;
		
		if(juve > 0.0)
		{
			float js <- juvenile_accumulator;
			if(juve >= js)
			{
				res_qte_juve <- kill_all(false, true,false)[1];
	
			}
			else
			{
				int i <-0;
				loop while:i<length(juvenile_day_acc)
				{
					float to_kill <- juvenile_day_acc[i]/js * juve;
					res_qte_juve <- res_qte_juve + to_kill;
					juvenile_day_acc[i] <- juvenile_day_acc[i] - to_kill;
					juvenile_day_acc[i] <- juvenile_day_acc[i]>0?juvenile_day_acc[i]:0;
					i <- i +1;
				}
			}
		}
		if(young > 0.0)
		{
			float  js<- younger_accumulator;
			if(young >=js)
			{
				res_qte_young <- kill_all(false, false,true)[2];
			}
			else
			{
				int i <-0;
				loop while:i<length(young_born_day_acc)
				{
					float to_kill <-  juvenile_day_acc[i]/js * young;
					res_qte_young <- res_qte_young + to_kill ;
					young_born_day_acc[i] <- young_born_day_acc[i] - to_kill;
					young_born_day_acc[i] <- young_born_day_acc[i]>0?young_born_day_acc[i]:0;
					i <- i +1;
				}
									
			}
		}
		list<float> to_ret<- convert_ind2weight(res_qte_adult,res_qte_juve,res_qte_young);
		return to_ret;
	}

	/**
	 * population management
	 */
	reflex transfertPopulations 
	{
		//Retrieve population that must be transfered to the upper age category
		float new_juvenile <- pull_data(young_born_day_acc,younger_compartment_max_age);
		float new_adult <- pull_data(juvenile_day_acc,juvenile_compartment_max_age);



		//Population transfert to the upper age category
		juvenile_day_acc <- juvenile_day_acc + new_juvenile;
		adult_accumulator <- adult_accumulator + new_adult;
		
		//Population growth
		float temp_repo <- current_date.month>= open_season and current_date.month<= close_season?  adult_grow_rate:0;
		float juvenile_sum <- sum(juvenile_day_acc);
			
		float next_step_younger_accumulator <- temp_repo*step * (adult_accumulator + juvenile_sum) *(1-total_population/carring_capacity);
		if (adult_accumulator !=0)
		{
//			add next_step_younger_accumulator to:young_born_day_acc;
			young_born_day_acc <- young_born_day_acc + next_step_younger_accumulator;
		}
		else
		{
			young_born_day_acc <- young_born_day_acc + next_step_younger_accumulator;
		}
		/**
		 * Substract adult, juvenile and younger dead voles from the population
		 */		
		float young_sum <- sum(young_born_day_acc);
		float to_young_death <- young_sum=0?0:young_sum*younger_death_rate*step/length(young_born_day_acc);
		float to_juvenile_death <-juvenile_sum=0?0: juvenile_sum*juvenile_death_rate*step/length(juvenile_day_acc);
		juvenile_day_acc <- juvenile_day_acc collect (each<to_juvenile_death?0: each - to_juvenile_death);
		young_born_day_acc <- young_born_day_acc collect (each<to_young_death?0:each - to_young_death);
		adult_accumulator <- adult_accumulator -  adult_accumulator*adult_death_rate*step;
	}
	
	/**
	 *	Create vole agent to initiate vole spatial spread  
	 */	
	reflex  initiate_vole_move when:  (carring_capacity != 0 and total_population/carring_capacity>dispersion_threshold)
	{
		//determine the number of agent to create
		int number_individual_move <-round(total_population-(carring_capacity*dispersion_threshold));
		number_individual_move <- juvenile_accumulator < number_individual_move ? int(juvenile_accumulator):number_individual_move; 
	
		float len <- juvenile_accumulator;
		
		//create vole agent.
		if(number_individual_move > 0)
		{
			juvenile_day_acc <- juvenile_day_acc collect (each - number_individual_move * each/len );
			point loc <- location;
		
			create voles number:number_individual_move
			{
				// select one vole in the juvenile age list
				juvenile_age_id <- rnd(length(myself.juvenile_day_acc)-1);
				// locate the agent at the current cell
				location <- loc;
				my_current_cell <- myself;
				current_place <- my_current_cell.my_mnt_cell; 
				//is the current vole agent a long disperser?
				crazy_voles <- flip(crazy_vole_rate);
			}
				
			
		}
	}
	
}
/**
 * starting place of the population spread
 */
species plots
{
	aspect default
	{
		draw shape color:#red;
	}
}

/* experiment to play */
experiment campagnols type: gui {
	/**
	 * parameters
	 */
	parameter "simnum" var: simnum <- 0 ;
	parameter "younger_compartment_change_rate" var: younger_compartment_max_age <- 18#day ;
	parameter "juvenile_compartment_change_rate" var:juvenile_compartment_max_age<- (8*7)#day ;
	parameter "adult_grow_rate" var:adult_grow_rate <- 0.0165/#day;
	parameter "younger_death_rate" var:younger_death_rate; // <- 0.003/#day;
	parameter "juvenile_death_rate" var:juvenile_death_rate; //<-0.003/#day;
	parameter "adult_death_rate" var:adult_death_rate; //<- 0.003/#day;
	parameter "dispersion_threshold" var:dispersion_threshold<- 0.2;
	parameter "default_vole_carring_capacity" var:default_vole_carrying_capacity<- 500/(100#m*100#m);
	parameter "crazy_vole_rate" var:crazy_vole_rate <- 0.1;
	parameter "open_season" var:open_season <- 4;
	parameter "close_season" var:close_season <- 10;
	parameter "min_distance_crazy_vole" var:min_distance_crazy_vole <- 100#m;
	parameter "max_distance_crazy_vole" var:max_distance_crazy_vole <- 3000#m;
	parameter "buffer_size" var:buffer_size  <- 500#m;
/**
 * outputs of the model
 */	
	output {
		monitor "discovered_area_year" value:discovered_area_year;
		monitor "densite_monitor" value:densite_monitor;
		monitor "densite_occuped_monitor" value:densite_occuped_monitor;
		monitor "population_young_monitor" value:population_young_monitor;
		monitor "population_juvenile_monitor" value:population_juvenile_monitor;
		monitor "population_adult_monitor" value:population_adult_monitor;
		monitor "max_cell_population" value:max_cell_population;
		monitor "long_distance_dispersing_voles" value:long_distance_dispersing_voles;
		monitor "short_distance_dispersing_voles" value:short_distance_dispersing_voles;
		monitor "current_date" value:""+current_date.year+"-"+current_date.month+"-"+current_date.day;
		monitor "discovered_area" value:discovered_area;

		display predating_map
		{
			species death_factor aspect:base transparency:0.7;
		}

		display discovering_map
		{
			grid cells;
			species zone_information aspect:base;
			species death_factor aspect:discovering_aspect transparency:0.5;
		}

		display campagnols_map 
		{
			species plots aspect:default;
			grid cells;
			species voles aspect:default;
			species voles_start_point aspect:default;

		}
		
		display campagnols_population
		{
			chart "Species evolution" type: series  {
				data "younger" value: sum(cells collect(each.younger_accumulator)) color: #blue ;
				data "juvenile" value: sum(cells collect(each.juvenile_accumulator)) color: #green ;
				data "adult" value: sum(cells collect(each.adult_accumulator)) color: #red ;
				data "cumul" value: (sum(cells collect(each.total_population)))  color: #pink ;

			}
		}

		display campagnols_population_percent
		{
			chart "Species evolution percent" type: series  {
				data "younger" value: pop_younger/pop_total * 100  color: #blue ;
				data "juvenile" value: pop_juvenile/pop_total * 100  color: #green ;
				data "adult" value: pop_adult/pop_total * 100  color: #red ;
			}
		}
	}
}
