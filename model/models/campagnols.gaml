/**
 *  campagnols
 *  Author: nicolas
 *  Description: 
 */

model campagnols

global {
	file space_grid_file <- file("../includes/DEM_Romanche_LIIres50.asc");
	file emprise <- file("../includes/emprise1.shp");
	file plot_shape <- file("../includes/Mailles_2014_reduit.shp");
	
	file starting_shape <- file("../includes/historique_1998.shp");
	
	//file starting_shape <- file("../includes/startingPoint_proj.shp");
	map colors <- map([1::rgb([178, 180, 176]), 2::rgb([246, 111, 0]), 3::rgb([107, 0, 0]), 4::rgb([249, 0, 255]), 5::rgb([144, 96, 22]), 6::rgb([255, 255, 86]), 7::rgb([19, 114, 38]), 8::rgb("black"), 9::rgb([107, 94, 255]), 10::rgb([43, 255, 255])]);
	
	geometry shape <- envelope(emprise);
	list<cells> simulated_plots <- [];
	
	
	int simnum<-0;
	// nombre d'adulte / nb juvénile
	
	float younger_compartment_max_age <- 18#day ;
	float juvenile_compartment_max_age <- (8*7)#day ;
	
	float adult_grow_rate <- 0.0165/#day; //0.38/#day;
	
	float younger_death_rate <- 0.01/#day;
	float juvenile_death_rate <- 0.01/#day;
	float adult_death_rate <- 0.02/#day;
	
	float dispersion_threshold <- 0.2;
	
	float default_vole_carrying_capacity <- 500/(100#m*100#m);
	
	float crazy_vole_rate <- 0.1;
	float min_distance_crazy_vole <- 3#km;
	float max_distance_crazy_vole <- 6#km;
	
	int open_season <- 4;
	int close_season <- 10;
	
	//int current_month <- 1 update: (time mod #year) div #month + 1;
	
	float born_accumulator <- 0.0;
	float death_accumulator <-0.0;
	
	float pop_total <- 1 update: (sum(cells collect(each.total_population)));
	float pop_adult <- 0 update: sum(cells collect(each.adult_accumulator));
	float pop_juvenile <- 0 update: sum(cells collect(each.juvenile_accumulator));
	float pop_younger <- 0 update:sum(cells collect(each.younger_accumulator));
	float colonisation_speed <- 0.0;
	
	float discovered_area <- 0;
	float discovered_area_year <- 0;
	
	
	float densite_monitor <- 0 update:mean(cells collect(((each.total_population/each.shape.area)/default_vole_carrying_capacity)*100));
	float densite_occuped_monitor <- 0 update:mean(cells where(each.is_occupated) collect(((each.total_population/each.shape.area)/default_vole_carrying_capacity)*100));
	float population_young_monitor <- 0 update:sum(cells collect(each.younger_accumulator));
	float population_juvenile_monitor <- 0 update:sum(cells collect(each.juvenile_accumulator));
	float population_adult_monitor <- 0 update:sum(cells collect(each.adult_accumulator));
	float max_cell_population <- 0 update: max(cells collect(each.total_population));
	int long_distance_dispersing_voles <- 0 update:(voles count(each.crazy_voles));
	int short_distance_dispersing_voles <- 0 update:(voles count(!each.crazy_voles));
	int dispersing_voles -> {short_distance_dispersing_voles + long_distance_dispersing_voles };

	
//	float min_cell_younger_monitor <- 0 update:cells count(each.younger_accumulator !=0 and each.younger_accumulator = min([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
//	float min_cell_juvenile_monitor <- 0 update:cells count(each.juvenile_accumulator !=0 and each.juvenile_accumulator = min([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
//	float min_cell_adult_monitor <- 0 update:cells count(each.adult_accumulator !=0 and each.adult_accumulator = min([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
//	float max_cell_younger_monitor <-0 update:cells count(each.younger_accumulator !=0 and each.younger_accumulator = max([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
//	float max_cell_juvenile_monitor<- 0 update:cells count(each.juvenile_accumulator !=0 and each.juvenile_accumulator = max([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
//	float max_cell_adult_monitor<- 0 update:cells count(each.adult_accumulator !=0 and each.adult_accumulator = max([each.younger_accumulator,each.juvenile_accumulator,each.adult_accumulator]));
	
	float coeff_go_down <- 2/3;
	//date starting_date <- date(5#month);
	
	date dstart <- date([1998, 4,1,0,0,0]);
	
	init
	{
		//shape <- envelope(plot_shape);
		step <- 1#day;
		starting_date <- dstart;
		//save world to:"../includes/projection.shp" type:"shp";
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
		
		
		point orig <- point({908115.991,2017912.707});
		create voles_start_point from: starting_shape //with:[younger_init::int(read("younger")),juvenile_init::int(read("juvenile")),adult_init::int(read("adult"))]
		{
			adult_init <- 100;
			younger_init <- 0;
			juvenile_init <- 0;
			
			write "juvenile " + juvenile_init;
			write "younger" + younger_init;
			write "grow " + adult_grow_rate;
		}
		
		ask voles_start_point
		{			
			ask cells overlapping(self)
			{
				//adult_accumulator <- 30.0;
				color <- #red; 
				adult_accumulator <- float(myself.adult_init);
				int i <- 0;
				young_born_day_acc<- [];				
			}
			
			//write "coucou" + younger_born_day at "cdsf";
			
		}
		
		ask mnt
		{
				associated_cell <- one_of(cells overlapping self);
				associated_cell.my_mnt_cell <- self;
				//neighbours_mnt_cell <- (self neighbours_at 1 ) where(each.height <=self.height );
				
				neighbours_mnt_cell <- (self neighbours_at 1 ) ;
				int down <- neighbours_mnt_cell count(each.height <=self.height);
				sum_weight <- coeff_go_down * down + (1-coeff_go_down) * (length(neighbours_mnt_cell) - down);
				
				
				//list <mnt> places <- (self neighbours_at 1 ) where(each.height <=self.height )
				
		}
		
		ask simulated_plots
		{
			my_far_cell <- simulated_plots overlapping (max_distance_crazy_vole around circle(min_distance_crazy_vole));
			
		}
		
		
	}
	
	reflex no_vol when: false and population_young_monitor + population_juvenile_monitor + population_adult_monitor <= 0
	{
		do halt;
	}
	
	reflex new_year when:  dstart.month = current_date.month and dstart.day = current_date.day and current_date.year != dstart.year
	{
		float area_cell <- one_of(cells).shape.area;
		float discovered_cell <- cells count (each.is_occupated) * area_cell;
		discovered_area_year <- discovered_cell - discovered_area;
		discovered_area <- discovered_cell;
	} 
	
	reflex newsemester when: current_date.day=1 and (  current_date.month = 6 or current_date.month = 12)
	{
	save simulated_plots with:[younger_accumulator::"younger",juvenile_accumulator::"juvenile",adult_accumulator::"adult"] to:"outs/campagnols_"+simnum+"_"+cycle+".shp" type:"shp";
	} 
	
}


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

species voles //skills:[moving]
{
	float age <- 0.0;
	int id <- 0;
	bool crazy_voles<-true;
	mnt current_place <- nil;
	
	cells my_current_cell;
	
	
	reflex gogo_crazy when: crazy_voles
	{
		my_current_cell <- one_of(my_current_cell.my_far_cell); //one_of(simulated_plots overlapping (6#km around circle(3#km)));
		location <- my_current_cell.location;
	}
	reflex gogo when: !crazy_voles // and false
	{
		//list<mnt> places <- current_place.neighbours_mnt_cell; //neighbours_at 1 where(each.height <=current_place.height);
		mnt plc <- choose_cell(current_place,rnd(1000)/1000,coeff_go_down,current_place.sum_weight);    //length(places)>0?one_of(places):current_place;
		current_place <- plc;
		location <- current_place.location;
		my_current_cell <- current_place.associated_cell;// one_of(cells overlapping current_place);
		
	}
	
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
	

	reflex killed when:flip(0.5)
	{
		do die;
	}
	
	reflex live when: my_current_cell!=nil and my_current_cell.carring_capacity!= 0 and my_current_cell.total_population/my_current_cell.carring_capacity<dispersion_threshold
	{
		put my_current_cell.juvenile_day_acc at id + 1 at: id in: my_current_cell.juvenile_day_acc;  
		do die;
	}
	
	aspect default
	{
		draw circle(20#m) color:#blue;
	}
}

species voles_day_pop
{
	int creation_step;
	int population;
}

grid mnt   file: space_grid_file schedules:[] {
	float height <- grid_value;
	cells associated_cell <- nil;
	
	list< mnt> neighbours_mnt_cell  <- nil;
	float sum_weight <- 0.0;
	init {
		//write grid_value;
		//color <- colors at int(grid_value);
	}
	
}
grid cells cell_width:100#m cell_height:100#m schedules:  simulated_plots
{
	list<cells> my_far_cell<-[];
	mnt my_mnt_cell;
	list<float> young_born_day_acc <-[];
	list<float> juvenile_day_acc <-[];
	float younger_accumulator -> {sum(young_born_day_acc)};
	float juvenile_accumulator -> {sum(juvenile_day_acc)};
	float adult_accumulator <- 0.0;
	float total_population -> {sum(young_born_day_acc) + sum(juvenile_day_acc) +adult_accumulator };
	float adult_rate -> {adult_accumulator = 0? 0:adult_accumulator / total_population};
	float carring_capacity <- 0.0;
	//rgb color <- #white update: carring_capacity=0?#white:(total_population=0?#green:rgb(255,0,40));
	rgb color <- #white update: carring_capacity=0?#white:(total_population!=0?#red:#green);
	bool is_occupated -> {total_population>0};
	

	aspect grille
	{
		if(carring_capacity != 0)
		{
			draw shape color:#white;			
		}
	}
	
	
	
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
	
	list<float>  killed_voles_rated(float adult, float juvenile, float younger)
	{
		list<float> result <- [0.0,0.0,0.0];
		float juv <- juvenile_accumulator;
		if(juv < juvenile)
		{
			juvenile <- juvenile - juv;
			juvenile_day_acc<- [];
			result[0]<- juvenile;
		}
		else
		{
			juvenile_day_acc <- juvenile_day_acc accumulate (each - juvenile)	;			
		}
		float you <- younger_accumulator;
		if(you < younger)
		{
			younger <- younger - you;
			young_born_day_acc<- [];
			result[1]<- younger;
		}
		else
		{
			young_born_day_acc <- young_born_day_acc accumulate (each - younger)	;			
		}
		
		if(adult_accumulator < adult)
		{
			adult <- adult - adult_accumulator;
			adult_accumulator <- 0.0;
			result[2]<-adult;
		}
		else
		{
			adult_accumulator <- adult_accumulator - adult;			
		}
		return result;
	}
	
	action killed_voles(float nbInd)
	{
		float totPop <- total_population;
		
		float xx <- nbInd / totPop;
		juvenile_day_acc <- juvenile_day_acc accumulate (each - each * xx)	;
		young_born_day_acc <- young_born_day_acc accumulate (each - each * xx)	;
		adult_accumulator <- adult_accumulator - adult_accumulator*xx;
		
	}
	
	reflex transfertPopulations
	{
		float new_juvenile <- pull_data(young_born_day_acc,younger_compartment_max_age);
		float new_adult <- pull_data(juvenile_day_acc,juvenile_compartment_max_age);
		
		//write 'test juv ' + new_juvenile + " " + new_adult + " "+ current_month;
		
		juvenile_day_acc <- juvenile_day_acc + new_juvenile;
		adult_accumulator <- adult_accumulator + new_adult;
		float temp_repo <- current_date.month>= open_season and current_date.month<= close_season?  adult_grow_rate:0;
		float juvenile_sum <- sum(juvenile_day_acc);
		
		float next_step_younger_accumulator <- temp_repo*step * (adult_accumulator + juvenile_sum) *(1-total_population/carring_capacity);
		if (adult_accumulator !=0)
		{
					young_born_day_acc <- young_born_day_acc + next_step_younger_accumulator;
						
		}
		float young_sum <- sum(young_born_day_acc);
		float to_young_death <- young_sum=0?0:young_sum*younger_death_rate*step/length(young_born_day_acc);
		float to_juvenile_death <-juvenile_sum=0?0: juvenile_sum*juvenile_death_rate*step/length(juvenile_day_acc);
		
		//juvenile_day_acc <- juvenile_day_acc accumulate (each<to_juvenile_death?0: each - to_juvenile_death); 
		juvenile_day_acc <- juvenile_day_acc collect (each<to_juvenile_death?0: each - to_juvenile_death);
		//young_born_day_acc <- young_born_day_acc accumulate (each<to_young_death?0: each - to_young_death); 
		young_born_day_acc <- young_born_day_acc collect (each<to_young_death?0:each - to_young_death);
		adult_accumulator <- adult_accumulator -  adult_accumulator*adult_death_rate*step;
		
	}
	
	
	reflex  initiate_move when:  (carring_capacity != 0 and total_population/carring_capacity>dispersion_threshold)
	{
			
		int number_individual_move <-round(total_population-(carring_capacity*dispersion_threshold));
		 number_individual_move <- juvenile_accumulator < number_individual_move ? int(juvenile_accumulator):number_individual_move; 
		//float nb_juvenile_accumulator <- 10;//juvenile_accumulator - number_individual_move;
		
		float len <- juvenile_accumulator;
		if(number_individual_move > 0)
		{
			juvenile_day_acc <- juvenile_day_acc collect (each - number_individual_move * each/len );
			point loc <- location;
		
			create voles number:number_individual_move
			{
				
				id <- rnd(length(myself.juvenile_day_acc)-1);
				location <- loc;
				my_current_cell <- myself;
				current_place <- my_current_cell.my_mnt_cell; //one_of(mnt where(each.associated_cell=my_current_cell) );
				crazy_voles <- flip(crazy_vole_rate);
			}
				
			
		}
	}
	
}

species plots
{
		
	aspect default
	{
		draw shape color:#red;
	}
}



experiment campagnols type: gui {
	/** Insert here the definition of the input and output of the model */
	parameter "simnum" var: simnum <- 0 ;
	parameter "younger_compartment_change_rate" var: younger_compartment_max_age <- 18#day ;
	parameter "juvenile_compartment_change_rate" var:juvenile_compartment_max_age<- (8*7)#day ;
	parameter "adult_grow_rate" var:adult_grow_rate <- 0.0165/#day;
	parameter "younger_death_rate" var:younger_death_rate <- 0.003/#day;
	parameter "juvenile_death_rate" var:juvenile_death_rate<-0.003/#day;
	parameter "adult_death_rate" var:adult_death_rate<- 0.003/#day;
	parameter "dispersion_threshold" var:dispersion_threshold<- 0.2;
	parameter "default_vole_carring_capacity" var:default_vole_carrying_capacity<- 500/(100#m*100#m);
	parameter "crazy_vole_rate" var:crazy_vole_rate <- 0.1;
	parameter "open_season" var:open_season <- 4;
	parameter "close_season" var:close_season <- 10;
	parameter "min_distance_crazy_vole" var:min_distance_crazy_vole <- 100#m;
	parameter "max_distance_crazy_vole" var:max_distance_crazy_vole <- 3000#m;
	
	
	
	
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
//		monitor "min_cell_younger_monitor" value:min_cell_younger_monitor;
//		monitor "min_cell_juvenile_monitor" value:min_cell_juvenile_monitor;
//		monitor "min_cell_adult_monitor" value:min_cell_adult_monitor;
//		monitor "max_cell_younger_monitor" value:max_cell_younger_monitor;
//		monitor "max_cell_juvenile_monitor" value:max_cell_juvenile_monitor;
//		monitor "max_cell_adult_monitor" value:max_cell_adult_monitor;
		monitor "current_date" value:""+current_date.year+"-"+current_date.month+"-"+current_date.day;
		monitor "discovered_area" value:discovered_area;
	
		display campagnols_map 
		{
			species plots aspect:default;
			grid cells;
			species voles aspect:default;
			species voles_start_point aspect:default;

		}
		
	/*	display campagnols_population
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
 */
		
		
	}
}
