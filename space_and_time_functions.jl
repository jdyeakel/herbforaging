# these functions handle simulation mechanics like generating arrays or transforming meta-parameters into their various use cases. 

function generate_loop_array(number_species, number_scenarios, number_strategies,
        configurations)
    # generates an appropriate array for holding intermediate data based on spp list and metaparameters
    total = number_species * number_scenarios  * 
    number_strategies *configurations
    toy_array = SharedArray{Int64}(total, 4);
    
    array_index=1

    for species in 1:number_species, scenario in 1:number_scenarios,
        target in 1:number_strategies, config in 1:configurations
    
        toy_array[array_index, 1] = species
        toy_array[array_index, 2] = scenario
        toy_array[array_index, 3] = target
        toy_array[array_index, 4] = config
        array_index += 1
    end
    return toy_array
end


function update_compartments(t, mouth, gut, fat, rates, costs, resource_gain;
        cropping=false, chewing=false, travelling=false, resource=0)
   # this version of update_compartments is not in use. Current function below
   # This function updates mouth,gut,fat levels according to events, behavior state, and time passage
   # takes: time, current mouth, gut, and fat levels, intake mass
    
    # cropping puts food in the mouth
    if cropping == true && chewing == false
        mouth += rates[1][resource]; 
    end

    # chewing takes food from the mouth to the gut
    if cropping==false && chewing == true
       if mouth >= rates[2] # 
           gut += rates[2] * t;
            mouth -= rates[2] * t;
        elseif mouth >0 && mouth < beta
            mouth -= mouth;
            gut += mouth;
        end
    end


    if cropping==true && chewing == true
        return print("you can't chew and crop")
    end

   if gut >0 # food moves from gut to fat whenever time passes
        gut -= min(gut, rates[3]*t);
        fat += min(gut, rates[3]*t*16.7); #conversion g -> kj
    end

    # fat depletes whenever time passes at a rate that depends on movement state (basal/field)
    if travelling==false
        fat -= min(fat, costs[1] * t); #
    elseif travelling==true
        fat -= min(fat, costs[2] * t); #
    end
    #elseif fat <= 0
        #return print("you have died")

    return mouth, gut, fat

end


function pool_prop(resource_pool, index)
    # helper function for determing average chewing rate depending on poportional contribution of 
    # grass/browse to content of mouth
    prop = resource_pool[index]/sum(resource_pool)
    return prop
end

function update_compartments2(t, mouth, gut, fat, rates, costs, resource_gain;
        cropping=false, chewing=false, travelling=false, resource=0)

   #  This function updates mouth,gut,fat levels according to events, behavior state, and time passage
   # takes: time, current mouth, gut, and fat levels, intake mass
    
    #cropping - food enters the mouth
    if cropping == true && chewing == false
        mouth[resource] += min(rates[1][resource]*t, 1.0);
    end

    #chewing - food moves from mouth to gut
    if cropping==false && chewing == true
        
        grass_chew_rate = pool_prop(mouth, 1) * rates[2]
        browse_chew_rate = pool_prop(mouth, 2) * rates[2]
        
       if sum(mouth) >= (rates[2]*t)
           gut[1] += grass_chew_rate * t;
           gut[2] += browse_chew_rate * t;
           mouth[1] -= grass_chew_rate * t;
           mouth[2] -= browse_chew_rate * t;
        elseif sum(mouth) >0 && sum(mouth) < rates[2]*t
           gut[1] += mouth[1];
           gut[2] += mouth[2];
           mouth[1] -= mouth[1];
           mouth[2] -= mouth[2];
        end
    end


    if cropping==true && chewing == true
        return print("you can't chew and crop")
    end

   if sum(gut) > 0.0 # food moves from gut to fat whenever time passes
 
        grass_digest_rate = pool_prop(gut, 1) * rates[3]
        browse_digest_rate = pool_prop(gut, 2) * rates[3]
        gut[1] -= min(gut[1], grass_digest_rate *t);
        gut[2] -= min(gut[2], browse_digest_rate *t);
        fat += (grass_digest_rate*resource_gain[1] + browse_digest_rate*resource_gain[2])*t
        
     
    end

    # fat depletes whenver time passes according to movement state (basal/field)
    if travelling==false
        fat -= min(fat, costs[1] * t); #
    elseif travelling==true
        fat -= min(fat, costs[2] * t); #
    end
    #elseif fat <= 0
        #return print("you have died")
    return mouth, gut, fat

end



function catch_food(chosen_resource, num_succ)
    # Do you catch the food?
    # will currently always catch food because bushes dont' run
    catch_food = rand();
    if catch_food < catch_success[chosen_resource]
    #You caught it! Pat yourself on the back

    # indexing here works, i think, b/c weight of target 1 is 0, so this never
    # comes up. If the weighting of target 1 changes, this will break
    num_succ[chosen_resource] += 1;
    end
        return num_succ
end


function forage(strategy_id, target_weight, target, resource_stats, bernouli_dist, gamma_dist)
    # helper function for "travel"
    # takes resource landscape and strategy, selects a resource to forage for
    let nearest_resource = 0,
    #t = 0.0, # start the clock for the day
    nearest_distance = 0.0 # init the distance of nearest resource
    distances = dist_to_resources(resource_stats, gamma_dist);
    nearest_distance, nearest_resource = find_nearest_resource(resource_stats, distances)
    chosen_resource = which_resource(strategy_id, target_weight,
            target, nearest_resource, bernouli_dist)
    chosen_distance = distances[chosen_resource]

        return chosen_resource, chosen_distance
    end
        #return chosen_resource

end

function travel(strategy_id, target_weight, target, resource_stats, velocity, t, t_max, bernouli_dist, gamma_dist)
    # selects next food item with "forage" and travels there
    resource, distance = forage(strategy_id, target_weight, target, resource_stats, bernouli_dist, gamma_dist)
    travel_time = distance/velocity;

    if t_max >= travel_time + t # check to make sure there is enough time left to get there
        return resource, travel_time
        elseif t_max < travel_time + t # if not, reroll
            travel(strategy_id, target_weight, target, resource_stats, velocity, t, t_max,  bernouli_dist, gamma_dist)
        else
            return print("something is awry in travel")
        end
    end


# checked basic function with res =  [[1 2]; [3 4]; [5 6]]
# haven't done full check-expect





function resource_selection(strategy_id,target_weight, target, resource_stats, bernouli_dist, gamma_dist)
    # pulls distance to resources, nearest resource, and which_res into 1 function call

    dist_to_resource = dist_to_resources(resource_stats, gamma_dist);
    nearest_distance, nearest_resource = find_nearest_resource(resource_stats, dist_to_resource);
            chosen_resource = which_resource(strategy_id, target_weight, target, nearest_resource, bernouli_dist);

    return chosen_resource, dist_to_resource[chosen_resource]
end


function which_resource(strategy_id, target_weight, target, nearest_resource, bernoulli_dist)
    # chooses between a targeted resource and the nearest resource
    # produces single Int value resource (1 or 2)
    # weighted coin flip, higher weight, more likely 1
    #bernoulli_dist = Bernoulli(1-target_weight[target]);
    draw = rand(bernoulli_dist);

    #If we draw a 0, the consumer will move towards the TARGETED RESOURCE
    if draw == 0 # if coinflip=0 (more likely if lower weighting)
    # which resource is targeted?
    chosen_resource = strategy_id[target];

    #If we draw a 1, the consumer will move towards the CLOSEST RESOURCE
        elseif draw == 1
        chosen_resource = nearest_resource
    end

    return chosen_resource #, draw, strategy_id[target]
end

function generate_bern_array(target_weight)
    # takes weightings for targeting decisions and returns array of bernoulli 
    # distrobutions representing those weights
    number_weights = length(target_weight);
    bern_array = []
    for i = 1:number_weights
    
        push!(bern_array, Bernoulli(1-target_weight[i]))
    end
    return bern_array
end


function generate_gamma_array(resource_scenarios)
    # takes the plant resource distribution scenarios and generates an array of gamma distributions
    # representing each scenario
    number_scenarios = length(resource_scenarios);
    gamma_array = []
    for i = 1:number_scenarios

        resource_stats = resource_scenarios[i]
        push!(gamma_array, resource_distributions(resource_stats[1,:], resource_stats[2,:]))
    end
    return gamma_array
end

function resource_distributions(observed_mean, observed_variance)
    # this function takes observed landscape mean and variance for plant resources and translates them 
    # into  alpha and c parameters for gamma distributions
    # developed in concert with Uttam
    c = observed_mean ./ (observed_variance .- observed_mean); # [float] [g/m2 / (g2/m4 * g/m2) = g/m2 / g3/m6 = m4/g2]
    alpha = c .* observed_mean;   # [float] [m4/g2 * g/m2 = m2/g]
    #lambda = alpha ./ c; # [float vector] [g/m2]
    r = alpha;
    p = 1.0 ./ (1.0 .+ c);
    #Define the negative binomial distribution
    #neg_bin = NegativeBinomial(r, p)

    #Define the gamma distribution
    gamma_dist = Gamma.(alpha, 1.0 ./ c); #mean = alpha * (1/c)

    return gamma_dist
end



function dist_to_resources(resource_stats, gamma_dist)
    # function for building list of distances to each resource
    # produces distances for [res1, res2]
    number_resources = length(resource_stats[:,1])
    #Draw distances to each resource (meters)
    #gamma_dist = resource_distributions(resource_stats[1,:], resource_stats[2,:]);
    distance_to_resource = zeros(Float64, number_resources);

    for i = 1:number_resources # for each resource
    # rand can be given a random number generator
    # this is the 2-step distance draw because you can't just build a neg bin
        distance_to_resource[i] = rand(Exponential(1.0 / rand(gamma_dist[i]))); # in meters
          end

    return distance_to_resource
end



function find_nearest_resource(resource_stats, dist_to_resource)
    # What resource is the closest?
    # produces tuple of (distance, resource_id)
    distance_tuple = findmin(dist_to_resource); # touple b/c (value, index)
    nearest_distance = distance_tuple[1]; # value of touple
    nearest_resource = distance_tuple[2]; # index of touple
    return nearest_distance, nearest_resource
end


function strategy_matrix_translation(strat_matrix)
    # takes a matrix of optimal strategies and returns one with plotting values for the plot with a 
    # matrix of green/yellow restults by landscape and bodysize
    dim1 = size(strat_matrix)[1];
    dim2 = size(strat_matrix)[2];
    translated_matrix = zeros(Int64, dim1, dim2); 
    
    for i in 1:dim1, j in 1:dim2
        if strat_matrix[i, j] == 1
            translated_matrix[i,j] = 5
            
            elseif strat_matrix[i, j] == 2.0
            translated_matrix[i,j] = 4
            
            elseif strat_matrix[i, j] == 3.0
            translated_matrix[i,j] = 3
            
            elseif strat_matrix[i, j] == 4.0
            translated_matrix[i,j] = 2
            
            elseif strat_matrix[i, j] == 5.0
            translated_matrix[i,j] = 1
            
            elseif strat_matrix[i, j] == 6.0
            translated_matrix[i,j] = 5
            
            elseif strat_matrix[i, j] == 7.0
            translated_matrix[i,j] = 6
            
            elseif strat_matrix[i, j] == 8.0
            translated_matrix[i,j] = 7
            
            elseif strat_matrix[i, j] == 9.0
            translated_matrix[i,j] = 8
            
            elseif strat_matrix[i, j] == 10.0
            translated_matrix[i,j] = 9
        end
    end
    
    return translated_matrix
    
    end    
    

#This is a note of the translation from simulation strategy to plotting strategy

#1,6 are random     5
#2 low pref graze   4
#3 med pref graze   3
#4 high pref graze  2
#5 graze only       1

#7 low prey browse  6
#8 med pref browse  7
#9 high pref browse 8
#10 browse only     9


function best_strategy(body_condition)
    # takes body condition data from the range of runs
    # returns which foraging strategy produced the best results in each scenario
    # try mean and std of configs 
    for i=1:number_scenarios, j=1:number_strategies
        # take some kind of average over the configs for each scene/strat combo
        # probably also check std
        # compare the different strategies for each scenario
        # identify best strategy
        # report which best strategy
        
        
    end
        
end