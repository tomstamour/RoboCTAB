# RoboCTAB DNA extraction protocol for the OT-2 robot from Opentrons. Read the full method article here: https://www.mdpi.com/2223-7747/14/15/2263
# Automated version of the CTAB DNA extraction method, adapted from J. J., and J. L. Doyle. (1987). 
# For the automated RoboCTAB implementation, please cite: "Recent advances in robotic DNA extraction have improved throughput for large-scale genotyping studies (Boucher St-Amour et al., 2025)."

############################################################
# Enter your values at lines 7 to 24 
first_column_plate_1 = 1                                #The first column in plate 1 for which you have samples.
last_column_plate_1  = 12                               #The last column in plate 1 for which you have samples.
first_column_plate_2 = 1
last_column_plate_2  = 12
first_column_plate_3 = 1
last_column_plate_3  = 12
first_column_plate_4 = 1
last_column_plate_4  = 12     

elution_buffer_volume = 40                              # Set the volume of Elution buffer (uL) that will be used to resuspend your DNA after isolation.
tipsbox = 'opentrons_96_tiprack_300ul'                  # Enter the API name of you tipsrack labware (300 uL without filters needed). 
samples_plate_type = '1.2ml_simport_vwr_t1102_96well'   # Enter the API name for the plate.s in which the samples were collected. This protocol is optimized for 1.2ml 96-well plates.
final_plate_type = '1.2ml_simport_vwr_t1102_96well'     # Enter the API name for the final plate.s in which the aquous phase (DNA) will be transfered too. This protocol is optimized for 1.2ml 96-well plates.
reservoir_type = 'agilent_1_reservoir_290ml'            # The reservoir must have a minimum capacity of 200ml.
chloroform_buffer_mixing = 'pipette_mixing'             # String: 'pipette_mixing' or 'bubble_mixing' or 'no_mixing'. When 'pipette_mixing' is selected, the robot will mix by doing inspirations and expulsions in the liquids. Select this option if your samples have 5mg or less of WELL GROUNDED material, otherwise the pipettes could block. When selecting 'bubble_mixing' the robot will blow air in the liquids. Select this method if your samples are over 5 mg and the grinding is poor. When selecting 'no_mixing', the robot will simply add the chloroform to the wells and emulsion/mixing will have to be done by either vortexing or inversions.
pipetteOff_isopropanol = False                          # If set to True, the robot will discard by pipetting-off the wash solution (isopropanol). If set to False the user will have to gently invert the plate to poor off and the discard the isopropanol.
pipetteOff_ethanol = False                              # If set to True, the robot will discard by pipetting-off the wash solution (ethanol). If set to False the user will have to gently invert the plate to poor off and the discard the ethanol.
distance_interstice_to_bottom = 14                      # To obtain this value, add 400uL of water to an empty tube of the "samples_plate_type" and measure the heigth (mm) of the liquid from the botom of the tube.
############################################################


# Running the Opentrons API
from opentrons import protocol_api
from opentrons import types
import math
import time

metadata = {'protocolName': 'RoboCTAB -- v1.1 --', 'apiLevel': '2.15'}

def get_values(*names):
    import json
    _all_values = json.loads("""{}""")
    return [_all_values[n] for n in names]

# Defining the total number of columns with samples in each plates
if first_column_plate_1 == 0 : number_of_columns_plate_1 = 0 
else: number_of_columns_plate_1 = (last_column_plate_1-first_column_plate_1 + 1)
if first_column_plate_2 == 0 : number_of_columns_plate_2 = 0 
else: number_of_columns_plate_2 = (last_column_plate_2-first_column_plate_2 + 1)
if first_column_plate_3 == 0 : number_of_columns_plate_3 = 0 
else: number_of_columns_plate_3 = (last_column_plate_3-first_column_plate_3 + 1)
if first_column_plate_4 == 0 : number_of_columns_plate_4 = 0 
else: number_of_columns_plate_4 = (last_column_plate_4-first_column_plate_4 + 1)
total_number_of_columns = number_of_columns_plate_1 + number_of_columns_plate_2 + number_of_columns_plate_3 + number_of_columns_plate_4

def run(ctx):
    
    # testing how many plates to extract
    if  first_column_plate_1 != 0 and \
        first_column_plate_2 == 0 and \
        first_column_plate_3 == 0 and \
        first_column_plate_4 == 0 : number_of_plates_to_extract = 1
    if  first_column_plate_1 != 0 and \
        first_column_plate_2 != 0 and \
        first_column_plate_3 == 0 and \
        first_column_plate_4 == 0 : number_of_plates_to_extract = 2
    if  first_column_plate_1 != 0 and \
        first_column_plate_2 != 0 and \
        first_column_plate_3 != 0 and \
        first_column_plate_4 == 0 : number_of_plates_to_extract = 3
    if  first_column_plate_1 != 0 and \
        first_column_plate_2 != 0 and \
        first_column_plate_3 != 0 and \
        first_column_plate_4 != 0 : number_of_plates_to_extract = 4


# Defining the functions executed in the protocols
    def water_distribute(all_samples_vector, volume, bottom_dispense_heigth, return_tip, source):

        nb_consecutive_dispense = 300 // volume
        nb_full_distribution = len(all_samples_vector) // nb_consecutive_dispense

        p300.default_speed = 400
        p300.pick_up_tip()

        # Defining the vectors for distribution
        count = 1
        while count <= nb_full_distribution + 1:

            if count == 1 :
                samples_start = 0
                samples_end = nb_consecutive_dispense
            
            elif count > 1 and count <= nb_full_distribution:
                 samples_start = (count * nb_consecutive_dispense) - nb_consecutive_dispense 
                 samples_end = samples_start + nb_consecutive_dispense 
            
            elif count == (nb_full_distribution + 1):
                samples_start = (nb_full_distribution * nb_consecutive_dispense)
                samples_end = len(all_samples_vector)     
            
            destinations = all_samples_vector[samples_start:samples_end]

            p300.default_speed = 400
            p300.aspirate(300, source.bottom(z = 2.5), rate = 1)
            
            dispense_count = 1
            for column in destinations:

                if dispense_count == 1:
                    volume_to_dispense = volume
                    p300.default_speed = 400
                else:
                    volume_to_dispense = volume + 10
                    p300.default_speed = 100
               
                p300.dispense(volume_to_dispense, location = column.bottom(z = bottom_dispense_heigth), rate = 1)
                
               
                p300.air_gap(10)

                dispense_count += 1

            p300.default_speed = 400
            p300.blow_out(location = source)
            
            count += 1


    def pre_grinding_water_dispense():
        p300.pick_up_tip()
        p300.default_speed = 400
        for plate in plates:
            for d in plate:
                p300.aspirate(50, water_reservoir_01.bottom(z = 2.5), rate = 1)
                p300.air_gap(5)
                p300.dispense(50 + 5, d.top(0), rate = 2)
                p300.blow_out(d.top(0))
        p300.return_tip()

    def ExtractionBuffer_dispense():
        p300.default_speed = 400
        #p300.pick_up_tip(tiprack_1.wells()[0])
        for plate in plates:
            for d in plate:
                # First dispense
                p300.aspirate(200, reservoir_01.bottom(z = 2), rate = 0.85)
                time.sleep(1.5)
                p300.air_gap(40)
                p300.dispense(240, d.top(2), rate = 2)
                p300.aspirate(20, d.top(2))
                p300.dispense(20, reservoir_01.top(z = 0))

                # Second dispense
                p300.aspirate(200, reservoir_01.bottom(z = 2), rate = 0.85)
                time.sleep(1.5)
                p300.air_gap(40)
                p300.dispense(240, d.top(2), rate = 2)
                p300.aspirate(20, d.top(2))
                p300.dispense(20, reservoir_01.top(z = 0))
        p300.return_tip()         
    

    def Chloroform_distribute(all_samples_vector, volume, top_dispense_heigth, return_tip, source):

        nb_consecutive_dispense = 300 // volume
        nb_full_distribution = len(all_samples_vector) // nb_consecutive_dispense
        
        p300.default_speed = 400
        p300.pick_up_tip()
        p300.aspirate(25, water_reservoir_01.bottom(z = 2.5), rate = 4) # Aspirating the water "Liquid-Cap"
        p300.aspirate(25, water_reservoir_01.top(z = 5), rate = 4)
        
        count = 1
        while count <= nb_full_distribution + 1:

            if count == 1 :
                samples_start = 0
                samples_end = nb_consecutive_dispense
            
            elif count > 1 and count <= nb_full_distribution:
                 samples_start = (count * nb_consecutive_dispense) - nb_consecutive_dispense 
                 samples_end = samples_start + nb_consecutive_dispense                     
            
            elif count == (nb_full_distribution + 1):
                samples_start = (nb_full_distribution * nb_consecutive_dispense)
                samples_end = len(all_samples_vector)     
            
            destinations = all_samples_vector[samples_start:samples_end]


            p300.default_speed = 200
            p300.aspirate(240, reservoir_01.bottom(z = 2.5), rate = 4)
            p300.aspirate(10, location = reservoir_01.top(z = 2.5))

            for column in destinations:

                p300.dispense(volume + 10, location = column.top(z = top_dispense_heigth), rate = 1)
                p300.air_gap(10)

            p300.default_speed = 400
            p300.dispense(volume = 10, location = source)
            
            count += 1

        if return_tip == True:
            p300.return_tip()
        else:
            p300.drop_tip()

    def dispensing_chloroform_and_pipetteMixing():
        
        # Loop through the plates
        for plate, tiprack in zip(plates, mixing_tipracks):
            
            # Loop through the wells in the plate
            for samples_wells, tips_column in zip(plate, tiprack):
                
                p300.default_speed = 400

                p300.pick_up_tip(location = tips_column)

                 # Aspirating the water "Liquid-Cap"
                p300.aspirate(40, water_reservoir_01.bottom(z = 2.5), rate = 4) 
                p300.aspirate(55, water_reservoir_01.top(z = 5), rate = 4) 

                # First dispensing
                p300.aspirate(200, reservoir_01.bottom(z = 2.5), rate = 2) 
                p300.aspirate(5, location = reservoir_01.top(z = 2.5)) 
                p300.dispense(195, location = samples_wells.top(z = 9), rate = 1)

                # Second dispensing
                p300.aspirate(193, reservoir_01.bottom(z = 2.5), rate = 2) 
                p300.aspirate(2, location = reservoir_01.top(z = 2.5))
                p300.dispense(195, location = samples_wells.bottom(z = distance_interstice_to_bottom + 6), rate = 1)
                
                # This loop will mix the chloroform and extraction buffer by blowing air in the liquids
                p300.aspirate(130, location = samples_wells.bottom(z = distance_interstice_to_bottom - 7), rate = 0.7)
                p300.default_speed = 400
                i = 0                                   
                while i < 5:  
                    p300.dispense(130, location = samples_wells.bottom(z = distance_interstice_to_bottom + 6), rate = 2)
                    p300.aspirate(130, location = samples_wells.bottom(z = distance_interstice_to_bottom - 7), rate = 0.7) 
                    i += 1

                p300.dispense(130, location = samples_wells.bottom(z = distance_interstice_to_bottom + 6), rate = 2)
                
                p300.touch_tip(samples_wells, v_offset = -3, radius = 1.2, speed = 40)

                p300.default_speed = 200
                p300.drop_tip(location = tips_column, home_after = False) # drop_tip() with no argument will drop the tips in the trash.
                p300.default_speed = 400

    time_dispensing_chloroform_and_bubbleMixing = 77             # Defining a time (s) variable to calculate the duration of the step


    def dispensing_chloroform_and_bubbleMixing():
        
        # Loop through the plates
        for plate, tiprack in zip(plates, mixing_tipracks):
            
            # Loop through the wells in the plate
            for samples_wells, tips_column in zip(plate, tiprack):
                
                p300.default_speed = 400

                p300.pick_up_tip(location = tips_column)

                 # Aspirating the water "Liquid-Cap"
                p300.aspirate(40, water_reservoir_01.bottom(z = 2.5), rate = 4) 
                p300.aspirate(55, water_reservoir_01.top(z = 5), rate = 4) 

                # First dispensing
                p300.aspirate(200, reservoir_01.bottom(z = 2.5), rate = 4) 
                p300.aspirate(5, location = reservoir_01.top(z = 2.5)) 
                p300.dispense(205, location = samples_wells.top(z = 9), rate = 1)

                # Second dispensing
                p300.aspirate(200, reservoir_01.bottom(z = 2.5), rate = 4) 
                p300.aspirate(5, location = reservoir_01.top(z = 2.5))
                p300.dispense(205, location = samples_wells.bottom(z = distance_interstice_to_bottom - 2), rate = 1)
                
                # This loop will mix the chloroform and extraction buffer by blowing air in the liquids
                p300.aspirate(205, location = samples_wells.top(z = 0), rate = 4)
                p300.default_speed = 400
                i = 0                                   
                while i < 10:  
                    p300.dispense(160, location = samples_wells.bottom(z = distance_interstice_to_bottom - 2), rate = 1)
                    p300.aspirate(160, location = samples_wells.top(z = 0), rate = 4) 
                    i += 1

                p300.dispense(180, location = samples_wells.bottom(z = distance_interstice_to_bottom - 2), rate = 1)
                
                p300.touch_tip(samples_wells, v_offset = -3, radius = 1.2, speed = 40)

                p300.default_speed = 200
                p300.drop_tip(location = tips_column, home_after = False) # drop_tip() with no argument will drop the tips in the trash.
                p300.default_speed = 400

    time_dispensing_chloroform_and_bubbleMixing = 77             # Defining a time (s) variable to calculate the duration of the step

    def dispensing_chloroform():
        p300.pick_up_tip()

        # Aspirating the water "Liquid-Cap"
        p300.aspirate(40, water_reservoir_01.bottom(z = 2.5), rate = 4) 
        p300.aspirate(55, water_reservoir_01.top(z = 5), rate = 4)
        
        # Loop through the plates
        for plate in plates:
            
            # Loop through the wells in the plate
            for samples_wells in plate:
            
                p300.default_speed = 400
                # First dispensing
                p300.aspirate(200, reservoir_01.bottom(z = 2.5), rate = 4) # Chloroform pipetting
                p300.aspirate(5, location = reservoir_01.top(z = 2.5)) # Air gap
                p300.dispense(205, location = samples_wells.top(z = 9), rate = 1)
                # Second dispensing
                p300.aspirate(200, reservoir_01.bottom(z = 2.5), rate = 4) # Chloroform pipetting
                p300.aspirate(5, location = reservoir_01.top(z = 2.5))
                p300.dispense(205, location = samples_wells.top(z = 9), rate = 1)
                
        p300.drop_tip()
    time_dispensing_chloroform = 7                      # Defining a time (s) variable to calculate the duration of the step

    volume_1 = 290 #
    volume_2 = 85 #
    
    def Supernatant_transfer(samples_plate, final_plate, transfer_tiprack):

        for source, destination, t in zip(samples_plate, final_plate, transfer_tiprack):
            
            p300.pick_up_tip(location = t)
            p300.move_to(source.top(-16), speed = 400)
            p300.move_to(source.bottom(z = distance_interstice_to_bottom + 3.25), speed = 15)                 # The z value is the distance to avoid touching the water-Chloroform interstice
            p300.aspirate(volume_1, source.bottom(z = distance_interstice_to_bottom + 3.25), rate = 0.6)      # The z value is the distance to avoid touching the water-Chloroform interstice
            p300.air_gap(10)

            p300.default_speed = 200
            p300.dispense(volume_1 + 10, destination.top(z = 1), rate = 1)            
            p300.blow_out(location = destination.top(z = -1))
            p300.touch_tip(location = destination, v_offset = -3, radius = 1.2, speed = 40)

            p300.move_to(source.top(-16), speed = 400)
            p300.move_to(source.bottom(z = distance_interstice_to_bottom + 2.5), speed = 7)                  # The z value is the distance to avoid touching the water-Chloroform interstice
            p300.aspirate(volume_2, source.bottom(z = distance_interstice_to_bottom + 2.5), rate = 0.2)    # The z value is the distance to avoid touching the water-Chloroform interstice
            p300.air_gap(10)

            p300.default_speed = 200
            p300.dispense(volume_2 + 10, destination.top(z = 1), rate = 2.5)
     
            p300.blow_out(location = destination.top(z = -1))
            p300.touch_tip(location =destination, v_offset = -3, radius = 1.2, speed = 40)
            p300.default_speed = 400
            p300.drop_tip(location = t, home_after = False)                            # This will return the tip in the tip rack at the same location were it was attached. This can be usefull if we would like to reuse the tips for removing isopropanol and ethanol from wash.

    time_supernatant_transfer = 45   

    def isopropanol_dispensing():
        isopropanol_volume = 295
        p300.pick_up_tip(tiprack_9["A1"])                  # transfer_tiprack_1 is always the firts tip rack for Isopropanol dispensing and supernanter stransfer
        p300.default_speed = 400
        
        for plate in final_plates:
            for d in plate:
                p300.aspirate(isopropanol_volume, reservoir_01.bottom(z = 2.5))
                p300.air_gap(5)
                p300.dispense(isopropanol_volume + 5, d.top(10))
                p300.aspirate(10, d.top(10))
                p300.dispense(volume = 10, location = reservoir_01.top(z = 4))
        p300.drop_tip()

    time_isopropanol_dispensing = 5                                 # Defining a time (s) variable to calculate the duration of the step

    def isopropanol_distribution(all_samples_vector, volume, top_dispense_heigth, return_tip, source):

        nb_consecutive_dispense = 300 // volume
        nb_full_distribution = len(all_samples_vector) // nb_consecutive_dispense
        
        p300.default_speed = 400
        p300.pick_up_tip(tiprack_9["A1"]) # Preventing a bug, it has to be define in this function
        count = 1
        while count <= nb_full_distribution + 1:

            if count == 1 :
                samples_start = 0
                samples_end = nb_consecutive_dispense
            
            elif count > 1 and count <= nb_full_distribution:
                 samples_start = (count * nb_consecutive_dispense) - nb_consecutive_dispense   
                 samples_end = samples_start + nb_consecutive_dispense                      
            
            elif count == (nb_full_distribution + 1):
                samples_start = (nb_full_distribution * nb_consecutive_dispense)
                samples_end = len(all_samples_vector)     
            
            destinations = all_samples_vector[samples_start:samples_end]

            p300.default_speed = 400
            p300.aspirate(290, source.bottom(z = 2.5), rate = 1)
            time.sleep(1)
            p300.air_gap(10)

            for column in destinations:

                p300.dispense(volume + 10, location = column.top(z = top_dispense_heigth), rate = 1)
                p300.air_gap(10)

            p300.default_speed = 400
            p300.blow_out(location = source)
            
            count += 1

        if return_tip == True:
            p300.return_tip()
        else:
            p300.drop_tip()

    def isopropanol_discarding(final_plate, transfer_tiprack):
        for s, t in zip(final_plate, transfer_tiprack):
            p300.default_speed = 400
            p300.pick_up_tip(location = t)

            # Doing a blow out on the side walls of the tiprak
            center_location = t.center()
            p300.blow_out(location = t.top(z = -8))
            p300.blow_out(center_location.move(types.Point(x = 1, y = 4, z = 21.5)))  
            p300.blow_out(center_location.move(types.Point(x = 1, y =-4, z = 21.5)))  
            p300.blow_out(center_location.move(types.Point(x = 4, y = 1, z = 21.5)))  
            p300.blow_out(center_location.move(types.Point(x =-4, y = 1, z = 21.5)))  
            p300.blow_out(center_location.move(types.Point(x = 0, y = 0, z = 21.5)))
              
            p300.aspirate(location =s.bottom(z = 8), volume = 295, rate = 0.7)
            p300.air_gap(5)
            p300.dispense(location = trash, volume = 300, rate = 2)
            p300.air_gap(5)


            p300.move_to(location = s.bottom(14))
            p300.dispense(location = s.bottom(14), volume = 5)
            p300.move_to(location = s.bottom(z = 4), speed = 30)            
            p300.aspirate(location =s.bottom(z = 4), volume = 295, rate = 0.4)
            p300.dispense(location = trash, volume = 300, rate = 2)
            p300.air_gap(5)

            p300.move_to(location = s.bottom(11))
            p300.dispense(location = s.bottom(11), volume = 5)
            p300.move_to(location = s.bottom(z = 1.75), speed = 30)            
            p300.aspirate(location =s.bottom(z = 1.75), volume = 70, rate = 0.1)
            p300.air_gap(5)

            p300.default_speed = 400
            p300.drop_tip(location = t, home_after = False)    # Drop_tip with no arguments will drop the tips in the trash.
    
    time_isopropanol_discarding = 55

    def ethanol_dispensing():
        p300.default_speed = 400
        p300.pick_up_tip(tiprack_9["A2"])
        for plate in final_plates:
            for f in plate:
                center_location = f.center()
                p300.aspirate(295, reservoir_01.bottom(z = 2.5), rate = 1.5)
                p300.air_gap(5)
                p300.dispense(300, center_location.move(types.Point(x = 1.25, y = 0, z = 22)), rate = 0.8) # Dispensing on the sidewall to avoid detachment of the DNA pellet at the bottom of the tubes.
                p300.blow_out(f.top(z = 4)) # Messy
                #p300.air_gap(5) # Useless
        p300.drop_tip()

    def ethanol_discarding(final_plate, transfer_tiprack):
        for s, t in zip(final_plate, transfer_tiprack):
            p300.pick_up_tip(location = t)
            p300.move_to(location = t.top(z = 0)) # Experimental to avoid the Homing after picking up the tip

            # Doing a blow out on the side walls of the tiprak
            center_location = t.center()
            p300.blow_out(location = t.top(z = -8))
            p300.blow_out(center_location.move(types.Point(x = 1, y = 4, z = 21)))    
            p300.blow_out(center_location.move(types.Point(x = 1, y =-4, z = 21)))    
            p300.blow_out(center_location.move(types.Point(x = 4, y = 1, z = 21)))    
            p300.blow_out(center_location.move(types.Point(x =-4, y = 1, z = 21)))    
            p300.blow_out(center_location.move(types.Point(x = 0, y = 0, z = 21))) 

            p300.aspirate(location =s.bottom(z = 5), volume = 200, rate = 0.7)
            
            p300.move_to(location = s.bottom(z = 3), speed = 20)
            p300.aspirate(location =s.bottom(z = 3), volume = 60, rate = 0.2)

            p300.move_to(location = s.bottom(z = 2), speed = 10)
            p300.aspirate(location =s.bottom(z = 2), volume = 35, rate = 0.1)
            p300.air_gap(5)
            
            p300.drop_tip(t, home_after = False)      # Drop_tip with no arguments will drop the tips in the trash.     


    time_ethanol_discarding = 32                       # Defining a time (s) variable to calculate the duration of the step

    def EBbuffer_dispensing():
            p300.default_speed = 400
            p300.pick_up_tip()    
            for plate in final_plates:   
                for f in plate:
                    p300.aspirate(elution_buffer_volume, reservoir_01.bottom(z = 1.5))
                    p300.air_gap(10)
                    p300.dispense(elution_buffer_volume + 10, f.top(0))
                    p300.blow_out(f.top(1))
            p300.drop_tip()

    def EB_distribute(all_samples_vector, volume, bottom_dispense_heigth, return_tip, source):

        nb_consecutive_dispense = 300 // volume
        nb_full_distribution = len(all_samples_vector) // nb_consecutive_dispense
        
        p300.default_speed = 400
        p300.pick_up_tip(tiprack_9["A3"])
        count = 1
        while count <= nb_full_distribution + 1:

            if count == 1 :
                samples_start = 0
                samples_end = nb_consecutive_dispense
            
            elif count > 1 and count <= nb_full_distribution:
                 samples_start = (count * nb_consecutive_dispense) - nb_consecutive_dispense 
                 samples_end = samples_start + nb_consecutive_dispense                       
            
            elif count == (nb_full_distribution + 1):
                samples_start = (nb_full_distribution * nb_consecutive_dispense)
                samples_end = len(all_samples_vector)     
            
            destinations = all_samples_vector[samples_start:samples_end]

            p300.default_speed = 400
            p300.aspirate(300, source.bottom(z = 1.75), rate = 1)
            
            dispense_count = 1
            for column in destinations:

                if dispense_count == 1:
                    volume_to_dispense = volume
                    p300.default_speed = 400
                else:
                    volume_to_dispense = volume + 10
                    p300.default_speed = 100
               
                p300.dispense(volume_to_dispense, location = column.bottom(z = bottom_dispense_heigth), rate = 1)
                p300.touch_tip(location = column, v_offset = 0.1, radius = 0.5, speed = 60) 
                p300.air_gap(10)

                dispense_count += 1

            p300.default_speed = 400
            p300.blow_out(location = source)
            
            count += 1

        if return_tip == True:
            p300.return_tip()
        else:
            p300.drop_tip()

    def truncate(n, decimals=0):                # This function is used to round decimal number for time calculation
        multiplier = 10**decimals
        return int(n * multiplier) / multiplier

#################################################################################################
#
#           REAGENTS                            REAGENTS                        REAGENTS 
#
#################################################################################################
    dead_volume_samples_number = 20
    samples_number = total_number_of_columns * 8
    
    Metabisulfite = round((samples_number + dead_volume_samples_number) * 0.0028, 3) # gram/sample
    PVPK29 = round((samples_number + dead_volume_samples_number) * 0.0056, 3) # gram/sample
    StockLysisSolution_A = round((samples_number + dead_volume_samples_number) * 0.2344, 1) # mL/sample
    StockLysisSolution_B = round((samples_number + dead_volume_samples_number)* 0.2344, 1) # mL/sample
    Sarkosyl  = round((samples_number + dead_volume_samples_number) * 0.0938,1) # mL/sample
    Rnase= int(round((samples_number + dead_volume_samples_number) * 0.2812, 0))  # uL/sample

    # Chloroform/Acool Isoamyl (24:1) calculator
    needed = (samples_number + dead_volume_samples_number) * 0.4 
    AlcoholIsoamyl = math.ceil(needed / 24)
    Chloroform = AlcoholIsoamyl * 24
    
    # Ethanol 70% calculator from 95% ethanol
    quantity_needed = int(round((samples_number + dead_volume_samples_number) * 0.3, 0))# in Ml
    C1 = 0.95  # C1 x V1 = C2 x V2
    C2 = 0.7
    V2 = quantity_needed
    V1 = int(round(C2 * V2 / C1, 1))
    volume_of_water1 = V2 - V1 # mL
    volume_of_ethanol95 = V1 # mL

#################################################################################################
#
#              COMMENTS                        COMMENTS                    COMMENTS
#
################################################################################################

    def comment_reagents_1(Metabisulfite, PVPK29, StockLysisSolution_A, StockLysisSolution_B, Sarkosyl, Rnase):
        return f'"Set water bath at 65°C. Prepare fresh working buffer : Stock Lysis Solution A {StockLysisSolution_A}ML + Stock Lysis Solution B {StockLysisSolution_B}ML + Sarkosyl {Sarkosyl}ML + Metabisulfite {Metabisulfite}g + PVP-K-29 {PVPK29}g. Store solution at 65C. Add {Rnase}uL of Rnase right before using solution'

    def comment_reagents_2(AlcoholIsoamyl, Chloroform, volume_of_water1, volume_of_ethanol95):
        return f'Chloroform-Isoamyl prep.: Chloroform {Chloroform}ML + Isoamyl Alcool {AlcoholIsoamyl}ML // Ethanol 70% prep.: DD Water {volume_of_water1}ML + Ethanol 95% {volume_of_ethanol95}ML (Store at -20°C)'

    comment_beginning = "Centrifuge the samples plates for 30s, then remove sealing tape and add metalic beads."
    
    def comment_1(p, sp, tr, wr):
        return f' Place samples plate {p} respectively on sites {sp} ---- Tipracks on sites {tr} ---- TE buffer reservoir on site {wr}'

    def comment_2(sp): 
        return f'START TE buffer  dispensing to samples on sites {sp}'
    
    comment_spinDown ="After grinding, centrifuge the plates for 30s then remove sealing tape."

    def comment_start_Chloro_dispensing(time_estimation):
        return f'START Chloroform dispensing to samples (~{time_estimation} minutes)'
    
    def comment_start_Supernatant_transfer(time_estimation):
        return f'START supernatant transfer (~{time_estimation} minutes)'
    
    def comment_start_Isopropanol_discarding(time_estimation):
        return f'START Isopropanol discarding (~{time_estimation} minutes)'
    
    def comment_start_Ethanol_discarding(time_estimation):
        return f'START Ethanol discarding (~{time_estimation} minutes)'
    
    def comment_decontaminate_chloroform(mixing_tipracks_sites):
        return f'Remove chloroform contaminated tipracks on sites {mixing_tipracks_sites} (dispose tips, clean tipracks)'
    
################################            ##    
#                              #          ####
# ONE plate (96-well) protocol #            ##
#                              #            ##   
################################           ####    
    if number_of_plates_to_extract == 1 :

        # Load Labware
        tiprack_1 = ctx.load_labware(tipsbox, '2')
        transfer_tiprack_1 = ctx.load_labware(tipsbox, '7')
        tiprack_9 = ctx.load_labware(tipsbox, '5')
        samples_plate_1 = ctx.load_labware(samples_plate_type, '1')
        final_plate_1 = ctx.load_labware(final_plate_type, '4')
        reservoir_1 = ctx.load_labware(reservoir_type, '3')
        water_reservoir = ctx.load_labware(reservoir_type, '6')
        trash = ctx.fixed_trash['A1']

        # Load Instrument
        p300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=[tiprack_1, transfer_tiprack_1, tiprack_9])
        p300.default_speed = 200

        # Subset only the columns with samples in your plates
        samples_P1 = samples_plate_1.rows()[0][first_column_plate_1 - 1:last_column_plate_1]
        plates = [samples_P1]
        all_samples = [well for plate in plates for well in plate]

        mixing_tiprack_1 = tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
        mixing_tipracks = [mixing_tiprack_1]

        final_plate_01 = final_plate_1.rows()[0][first_column_plate_1 - 1:last_column_plate_1]
        final_plates = [final_plate_01]
        all_final_plates = [well for plate in final_plates for well in plate]

        transfer_tiprack_01 = transfer_tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
        transfer_tipracks = [transfer_tiprack_01]

        reservoir_01 = reservoir_1.wells()[0]
        water_reservoir_01 = water_reservoir.wells()[0]


        #############################################################################################################
        #                                           DNA extraction actions
        #############################################################################################################
        
        samples_sites = '1'
        tipracks_sites = '2, 5 & 7'
        mixing_tipracks_sites = '2'
        water_reservoir_site ='6'
        reservoirs_sites = '3'

        ctx.pause(comment_reagents_1(Metabisulfite, PVPK29, StockLysisSolution_A, StockLysisSolution_B, Sarkosyl, Rnase))
        ctx.pause(comment_reagents_2(AlcoholIsoamyl, Chloroform, volume_of_water1, volume_of_ethanol95))
        ctx.pause('''Place full Isopropanol reservoir at -20°C''')
        ctx.pause(comment_beginning)
        ctx.pause(comment_1(samples_sites, samples_sites, tipracks_sites, water_reservoir_site))
        ctx.pause(comment_2(samples_sites))
        
        water_distribute(all_samples, volume = 50, bottom_dispense_heigth = 40, return_tip = True, source =  water_reservoir_01)

        ctx.pause('''Grind samples on a tyssus-lyser machine''')
        ctx.pause(comment_spinDown)
        ctx.pause('''Place samples plate back to site 1. Add Rnase to Extraction buffer and place in reservoir on site 3''')
        ctx.pause('''START post-grinding Extraction buffer dispensing to samples on site 1''')

        ExtractionBuffer_dispense()
        
        ctx.pause('''Seal plates with sealing tape and invert plates 10 times. Spin plates then remove sealing tape and incubate the plates (65C, 60 min)''')
        ctx.pause('''After incubation, place the plate back on site 1 and place Chloroform reservoir on site 3 and water reservoir on site 6''')
        time_estimation = str(truncate((time_dispensing_chloroform_and_bubbleMixing * total_number_of_columns / 60), 1))
        ctx.pause(comment_start_Chloro_dispensing(time_estimation))

        if chloroform_buffer_mixing == 'pipette_mixing':
            dispensing_chloroform_and_pipetteMixing()
            ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
            ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
        if chloroform_buffer_mixing == 'bubble_mixing':
            dispensing_chloroform_and_bubbleMixing()
            ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
            ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
        if chloroform_buffer_mixing == 'no_mixing':
            dispensing_chloroform()
            ctx.pause('''Mix (vortex carefully) then centrifugate the plate (6000rpm, 10 min).''')
        
        ctx.pause('''When centrifugation is done place the samples plate back on site 1 and place an empty plate (1.0ml 96-Deep well) on site 4 (label the plate)''')
        time_estimation = str(truncate(time_supernatant_transfer * total_number_of_columns / 60, 1))
        ctx.pause(comment_start_Supernatant_transfer(time_estimation))

        for samples, final, tiprack in zip(plates, final_plates, transfer_tipracks):
                    Supernatant_transfer(samples, final, tiprack)

        ctx.pause('''Remove Chloroform reservoir on site 3 and place cold Isopropanol reservoir on site 3''')
        ctx.pause('''START Isopropanol dispensing to plate on site 4''')

        isopropanol_dispensing()   

        ctx.pause('''At this point the plate (site 4) can be sealed and stored at -20C and the DNA extraction completed the next day (or later)'''  )
        ctx.pause('''Seal and invert 20x the plate on site 4 then centrifugate (10min, 6000 rpm)''') # and prepare cold 70 percent ethanol


        # Evaluating if the Isopropanol discarding is done by plate inversion or py pipetting off
        if pipetteOff_isopropanol == True:
            ctx.pause('''When centrifugation is done, place the plate back to site 4''')
            time_estimation = str(truncate(time_isopropanol_discarding * total_number_of_columns / 60, 1))
            ctx.pause(comment_start_Isopropanol_discarding(time_estimation))
        
            for final, tiprack in zip(final_plates, transfer_tipracks):
                isopropanol_discarding(final, tiprack)
        
        if pipetteOff_isopropanol == False:
            ctx.pause('''Gently invert the plate to poor off the supernatant then centrifuge the plate 10s''')
            ctx.pause('''Place plate back to site 4''')

        ctx.pause('''Remove Isopropanol reservoir on site 3, add Ethanol 70 percent reservoir to site 3''')
        ctx.pause('''START Ethanol dispensing to plate on site 4''')

        p300.starting_tip = tiprack_9.wells_by_name()['A1']

        ethanol_dispensing()

        ctx.pause('''Centrifugate the plate for 10 minutes at 6000 rpm''')

        # Ethanol discarding method (by inversion or py pipetting off)       
        if pipetteOff_ethanol == True:
            ctx.pause('''When centrifugation is done, place the plate back on site 4''')
            time_estimation = str(truncate(time_ethanol_discarding * total_number_of_columns / 60, 1))
            ctx.pause(comment_start_Ethanol_discarding(time_estimation))
            for final, tiprack in zip(final_plates, transfer_tipracks):
                ethanol_discarding(final, tiprack)
        if pipetteOff_ethanol == False:
            ctx.pause('''Gently invert the plate to poor off the supernatant then centrifuge the plate (6000rpm, 10sec)''')
        
        ctx.pause('''Evaporate ethanol (20min, 45C). Remove and clean tiprack on site 7 (place tips in trash)''')
        ctx.pause('''Prepare for Elution buffer dispensing: remove Ethanol reservoir on site 3 and place Elution buffer on site 3''')
        ctx.pause('''When evaporation is done, place plate back to site 4''')
        ctx.pause('''START Elution buffer dispensing''')

        p300.starting_tip = tiprack_9.wells_by_name()['A2']
        #EBbuffer_dispensing()
        EB_distribute(all_final_plates, volume = elution_buffer_volume, bottom_dispense_heigth = 16, return_tip = False, source = reservoir_01)

        ctx.pause('''Centrifuge plates at 4000rpm 5s. Seal plates with tape and store plates at 4C. DNA extraction completed''')
        ctx.comment('\n~~~~~~~~~~~~~~Protocol Complete~~~~~~~~~~~~~~\n')

    else:

###############################          ####
#                             #        ##    ##
# TWO plates 96-well protocol #             ##
#                             #           ##        
###############################         #######

        if number_of_plates_to_extract == 2 :
            
            # Load Labware
            tiprack_1 = ctx.load_labware(tipsbox, '6')
            tiprack_2 = ctx.load_labware(tipsbox, '9')
            transfer_tiprack_1 = ctx.load_labware(tipsbox, '7')
            transfer_tiprack_2 = ctx.load_labware(tipsbox, '8')
            tiprack_9 = ctx.load_labware(tipsbox, '10')
            samples_plate_1 = ctx.load_labware(samples_plate_type, '1')
            samples_plate_2 = ctx.load_labware(samples_plate_type, '2')
            final_plate_1 = ctx.load_labware(final_plate_type, '4')
            final_plate_2 = ctx.load_labware(final_plate_type, '5')
            reservoir_1 = ctx.load_labware(reservoir_type, '3')
            water_reservoir = ctx.load_labware(reservoir_type, '11')
            trash = ctx.fixed_trash['A1']

            # Load Instrument
            p300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=[tiprack_1, tiprack_2, transfer_tiprack_1, transfer_tiprack_2, tiprack_9])
            p300.default_speed = 200

            # Subset only the columns with samples in your plates
            samples_P1 = samples_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
            samples_P2 = samples_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
            plates = [samples_P1, samples_P2]
            all_samples = [well for plate in plates for well in plate]

            mixing_tiprack_1 = tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
            mixing_tiprack_2 = tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
            mixing_tipracks = [mixing_tiprack_1, mixing_tiprack_2]

            final_plate_01 = final_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
            final_plate_02 = final_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
            final_plates = [final_plate_01, final_plate_02]
            all_final_plates = [well for plate in final_plates for well in plate]

            transfer_tiprack_01 = transfer_tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
            transfer_tiprack_02 = transfer_tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
            transfer_tipracks = [transfer_tiprack_01, transfer_tiprack_02]
            
            reservoir_01 = reservoir_1.wells()[0]
            water_reservoir_01 = water_reservoir.wells()[0]

            #############################################################################################################
            #                                           DNA extraction actions
            #############################################################################################################
            
            samples_sites = '1 & 2'
            tipracks_sites = '6, 7, 8, 9 & 10'
            mixing_tipracks_sites = '6 & 9'
            water_reservoir_site ='11'
            reservoirs_sites = '3'

            ctx.pause(comment_reagents_1(Metabisulfite, PVPK29, StockLysisSolution_A, StockLysisSolution_B, Sarkosyl, Rnase))
            ctx.pause(comment_reagents_2(AlcoholIsoamyl, Chloroform, volume_of_water1, volume_of_ethanol95))
            ctx.pause('''Place full Isopropanol reservoir at -20°C''')
            ctx.pause(comment_beginning)
            ctx.pause(comment_1(samples_sites, samples_sites, tipracks_sites, water_reservoir_site))          
            ctx.pause(comment_2(samples_sites))
            
            water_distribute(all_samples, volume = 50, bottom_dispense_heigth = 40, return_tip = True, source =  water_reservoir_01)

            ctx.pause('''Grind samples on a tyssus-lyser machine''')
            ctx.pause(comment_spinDown)
            ctx.pause('''Place samples plates back to respective sites 1 & 2. Add Rnase to Extraction buffer and place in reservoir on site 3''')
            ctx.pause('''START dispensing of Extraction buffer to samples''')

            ExtractionBuffer_dispense()

            ctx.pause('''Seal plates with sealing tape and invert plates 10 times. Spin plates then remove sealing tape and incubate the plates (65C, 60 min)''')
            ctx.pause('''After incubation place the plates back to respective sites 1 & 2 and place Chloroform reservoir on site 3 and place water reservoir on site 11''')
            time_estimation = str(truncate((time_dispensing_chloroform_and_bubbleMixing * total_number_of_columns / 60), 1))
            ctx.pause(comment_start_Chloro_dispensing(time_estimation))

            if chloroform_buffer_mixing == 'pipette_mixing':
                dispensing_chloroform_and_pipetteMixing()
                ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
            if chloroform_buffer_mixing == 'bubble_mixing':
                dispensing_chloroform_and_bubbleMixing()
                ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
            if chloroform_buffer_mixing == 'no_mixing':
                dispensing_chloroform()
                ctx.pause('''Mix (vortex carefully) then centrifugate the plate (6000rpm, 10 min).''')

            ctx.pause('''When the 10 minutes centrifugation is done place the sample plates back to respective sites 1 & 2 and place empty plates (1.0ml 96-Deep well) on sites 4 and 5 (label the plates)''')
            
            time_estimation = str(truncate(time_supernatant_transfer * total_number_of_columns / 60, 1))
            ctx.pause(comment_start_Supernatant_transfer(time_estimation))
               
            for samples, final, tiprack in zip(plates, final_plates, transfer_tipracks):
                    Supernatant_transfer(samples, final, tiprack)

            ctx.pause('''Remove Chloroform reservoir on site 3 and place Isopropanol reservoir on site 3''')
            ctx.pause('''START Isopropanol dispensing to plates on site 4 & 5''')

            isopropanol_dispensing()

            ctx.pause('''At this point the plate (site 4 & 5) can be sealed and stored at -20C and the DNA extraction completed the next day (or later)'''  )
            ctx.pause('''Seal and invert 20x the plate on site 4 & 5 then centrifugate (10min, 6000 rpm)''') # and prepare cold 70 percent ethanol


            # Evaluating if the Isopropanol discarding is done by plate inversion or py pipetting off
            if pipetteOff_isopropanol == True:
                ctx.pause('''When centrifugation is done, place the plate back to respective sites 4 & 5''')
                time_estimation = str(truncate(time_isopropanol_discarding * total_number_of_columns / 60, 1))
                ctx.pause(comment_start_Isopropanol_discarding(time_estimation))

                for final, tiprack in zip(final_plates, transfer_tipracks):
                    isopropanol_discarding(final, tiprack)

            if pipetteOff_isopropanol == False:
                ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates 10s''')
                ctx.pause('''Place plates back to respective sites 4 & 5''')


            ctx.pause('''Remove Isopropanol reservoir on site 3, add Ethanol 70 percent reservoir to site 3''')
            ctx.pause('''START Ethanol dispensing to plates on site 4 & 5''')

            p300.starting_tip = tiprack_9.wells_by_name()['A1']

            ethanol_dispensing()

            ctx.pause('''Centrifugate the plates for 10 minutes at 6000 rpm''')

            # Ethanol discarding method (by inversion or py pipetting off)       
            if pipetteOff_ethanol == True:
                ctx.pause('''When centrifugation is done, place the plates back to respective sites 4 & 5''')
                time_estimation = str(truncate(time_ethanol_discarding * total_number_of_columns / 60, 1))
                ctx.pause(comment_start_Ethanol_discarding(time_estimation))
                for final, tiprack in zip(final_plates, transfer_tipracks):
                    ethanol_discarding(final, tiprack)
            if pipetteOff_ethanol == False:
                ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates (6000rpm, 10sec)''')

            ctx.pause('''Evaporate ethanol (20min, 45C). Remove and clean tipracks on sites 7 & 8 (place tips in trash)''')
            ctx.pause('''Prepare for Elution buffer dispensing: remove Ethanol reservoir on site 3 and place Elution buffer on site 3''')
            ctx.pause('''When evaporation is done, place plates back to respective sites 4 & 5''')
            ctx.pause('''START Elution buffer dispensing''')
            
            p300.starting_tip = tiprack_9.wells_by_name()['A2']
            #EBbuffer_dispensing()
            EB_distribute(all_final_plates, volume = elution_buffer_volume, bottom_dispense_heigth = 16, return_tip = False, source = reservoir_01)
            ctx.pause('''Centrifuge plates at 4000rpm 5s. Seal plates with tape and store plates at 4C. DNA extraction completed''')
            ctx.comment('\n~~~~~~~~~~~~~~Protocol Complete~~~~~~~~~~~~~~\n')
        
        else:
             
#################################            ###
#                               #           #  ##
# THREE plates 96-well protocol #             ##
#                               #           #  ##
#################################            ###
             
            if number_of_plates_to_extract == 3 :

                # Load Labware
                tiprack_1 = ctx.load_labware(tipsbox, '7')
                tiprack_2 = ctx.load_labware(tipsbox, '8')
                tiprack_3 = ctx.load_labware(tipsbox, '9')
                transfer_tiprack_1 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                transfer_tiprack_2 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                transfer_tiprack_3 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                tiprack_9 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                samples_plate_1 = ctx.load_labware(samples_plate_type, '1')
                samples_plate_2 = ctx.load_labware(samples_plate_type, '2')
                samples_plate_3 = ctx.load_labware(samples_plate_type, '3')
                final_plate_1 = ctx.load_labware(final_plate_type, '4')
                final_plate_2 = ctx.load_labware(final_plate_type, '5')
                final_plate_3 = ctx.load_labware(final_plate_type, '6')
                reservoir_1 = ctx.load_labware(reservoir_type, '10')
                water_reservoir = ctx.load_labware(reservoir_type, '11')
                trash = ctx.fixed_trash['A1']

                # Load Instrument
                p300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=[tiprack_1, tiprack_2, tiprack_3, transfer_tiprack_1, transfer_tiprack_2, transfer_tiprack_3, tiprack_9])
                p300.default_speed = 200

                # Subset only the columns with samples in your plates
                samples_P1 = samples_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                samples_P2 = samples_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                samples_P3 = samples_plate_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                plates = [samples_P1, samples_P2, samples_P3]                
                all_samples = [well for plate in plates for well in plate]

                mixing_tiprack_1 = tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                mixing_tiprack_2 = tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                mixing_tiprack_3 = tiprack_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                mixing_tipracks = [mixing_tiprack_1, mixing_tiprack_2, mixing_tiprack_3]

                final_plate_01 = final_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                final_plate_02 = final_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                final_plate_03 = final_plate_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                final_plates = [final_plate_01, final_plate_02, final_plate_03]
                all_final_plates = [well for plate in final_plates for well in plate]

                transfer_tiprack_01 = transfer_tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                transfer_tiprack_02 = transfer_tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                transfer_tiprack_03 = transfer_tiprack_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                transfer_tipracks = [transfer_tiprack_01, transfer_tiprack_02, transfer_tiprack_03]               
                
                reservoir_01 = reservoir_1.wells()[0]
                water_reservoir_01 = water_reservoir.wells()[0]

                #############################################################################################################
                #                                           DNA extraction actions
                #############################################################################################################
                
                samples_sites = '1, 2 & 3'
                tipracks_sites = '6, 7, 8, 9 & 10'
                mixing_tipracks_sites = '7, 8 & 9'
                water_reservoir_site ='11'
                reservoirs_sites = '3'

                ctx.pause(comment_reagents_1(Metabisulfite, PVPK29, StockLysisSolution_A, StockLysisSolution_B, Sarkosyl, Rnase))
                ctx.pause(comment_reagents_2(AlcoholIsoamyl, Chloroform, volume_of_water1, volume_of_ethanol95))
                ctx.pause('''Place full Isopropanol reservoir at -20°C''')
                ctx.pause(comment_beginning)
                
                ctx.pause('''Place Samples plates '1, 2 & 3' respectively on sites 1, 2 and 3''')
                ctx.pause('''Place tip racks on site 7, 8, 9 and place TE buffer reservoir on site 11''')
                ctx.pause('''START TE buffer  dispensing to samples''')
                
                water_distribute(all_samples, volume = 50, bottom_dispense_heigth = 40, return_tip = True, source =  water_reservoir_01)
                
                ctx.pause('''Grind samples on a tyssus-lyser machine''')
                ctx.pause(comment_spinDown)
                ctx.pause('''Place samples plates back to respective sites 1, 2 & 3. Add Rnase to Extraction buffer and place in reservoir on site 10''')
                ctx.pause('''START Extraction buffer dispenssing to samples on sites 1, 2 and 3''')
                
                ExtractionBuffer_dispense() 
                
                ctx.pause('''Seal plates with sealing tape and invert plates 10 times. Spin plates then remove sealing tape and incubate plates (65C, 60 min)''')
                ctx.pause('''After incubation, place the plates back to respective sites 1, 2, 3 and place Chloroform reservoir on site 10 and place water reservoir on site 11)''')
                time_estimation = str(truncate((time_dispensing_chloroform_and_bubbleMixing * total_number_of_columns / 60), 1))
                ctx.pause(comment_start_Chloro_dispensing(time_estimation))
                
                if chloroform_buffer_mixing == 'pipette_mixing':
                    dispensing_chloroform_and_pipetteMixing()
                    ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                    ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
                if chloroform_buffer_mixing == 'bubble_mixing':
                    dispensing_chloroform_and_bubbleMixing()
                    ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                    ctx.pause(comment_decontaminate_chloroform(mixing_tipracks_sites))
                if chloroform_buffer_mixing == 'no_mixing':
                    dispensing_chloroform()
                    ctx.pause('''Mix (vortex carefully) then centrifugate the plate (6000rpm, 10 min).''')
                
                for tiprack in [tiprack_1, tiprack_2, tiprack_3]:
                    ctx.move_labware(labware = tiprack, new_location = protocol_api.OFF_DECK)

                ctx.move_labware(labware = transfer_tiprack_1, new_location = 7)
                ctx.move_labware(labware = transfer_tiprack_2, new_location = 8)
                ctx.move_labware(labware = transfer_tiprack_3, new_location = 9)

                ctx.pause('''When the 10 minutes centrifugation is done place the sample plates back to respective sites 1, 2 & 3 and place empty plate (1.0ml 96-Deep well)s on sites 4, 5 & 6 (label the plates)''')
        
                time_estimation = str(truncate(time_supernatant_transfer * total_number_of_columns / 60, 1))
                ctx.pause(comment_start_Supernatant_transfer(time_estimation))

                for samples, final, tiprack in zip(plates, final_plates, transfer_tipracks):
                    Supernatant_transfer(samples, final, tiprack)

                ctx.pause('''Remove Chloroform reservoir on site 10 and place cold Isopropanol reservoir on site 10''')
                
                ctx.pause('''Remove reservoir on site 11 and place new tiprack on site 11''')
                ctx.move_labware(labware = water_reservoir, new_location = protocol_api.OFF_DECK)
                ctx.move_labware(labware = tiprack_9, new_location = 11)
                ctx.pause('''START Isopropanol dispensing to plates on site 4, 5 & 6''')

                isopropanol_dispensing()

                ctx.pause('''At this point the plate (sites 4, 5 & 6) can be sealed and stored at -20C and the DNA extraction completed the next day (or later)'''  )
                ctx.pause('''Seal and invert 20x the plate on sites 4, 5 & 6 then centrifugate (10min, 6000 rpm)''') # and prepare cold 70 percent ethanol

                # Evaluating if the Isopropanol discarding is done by plate inversion or py pipetting off
                if pipetteOff_isopropanol == True:
                    ctx.pause('''When centrifugation is done, place the plate back to respective sites 4, 5 & 6''')
                    time_estimation = str(truncate(time_isopropanol_discarding * total_number_of_columns / 60, 1))
                    ctx.pause(comment_start_Isopropanol_discarding(time_estimation))

                    for final, tiprack in zip(final_plates, transfer_tipracks):
                        isopropanol_discarding(final, tiprack)

                if pipetteOff_isopropanol == False:
                    ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates 10s''')
                    ctx.pause('''Place plates back to respective sites 4, 5 & 6''')

                ctx.pause('''Remove Isopropanol reservoir on site 10, add Ethanol 70 percent reservoir to site 10''')
                ctx.pause('''START Ethanol dispensing to plates on site 4, 5, 6''')

                ethanol_dispensing()

                ctx.pause('''Centrifugate the plates for 10 minutes at 6000 rpm''')

                # Ethanol discarding method (by inversion or py pipetting off)       
                if pipetteOff_ethanol == True:
                    ctx.pause('''When centrifugation is done, place the plates back to respective sites 4, 5 & 6''')
                    time_estimation = str(truncate(time_ethanol_discarding * total_number_of_columns / 60, 1))
                    ctx.pause(comment_start_Ethanol_discarding(time_estimation))
                    for final, tiprack in zip(final_plates, transfer_tipracks):
                        ethanol_discarding(final, tiprack)
                if pipetteOff_ethanol == False:
                    ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates (6000rpm, 10sec)''')

                ctx.pause('''Evaporate ethanol (20min, 45C). Remove and clean tipracks on sites 7, 8 & 9 (place tips in trash)''')
                ctx.pause('''Prepare for Elution buffer dispensing: remove Ethanol reservoir on site 10 and place Elution buffer on site 10''')
                ctx.pause('''When evaporation is done, place plates back to respective sites 4, 5 & 6''')
                ctx.pause('''START Elution buffer dispensing''')

                EB_distribute(all_final_plates, volume = elution_buffer_volume, bottom_dispense_heigth = 16, return_tip = False, source = reservoir_01)
                ctx.pause('''Centrifuge plates at 4000rpm 5s. Seal plates with tape and store plates at 4C. DNA extraction completed''')

            else:
                
#################################             ##    ##
#                               #             ##    ##
# FOUR plates 96-well protocol  #             ########      
#                               #                   ##
#################################                   ##
                 
                 if number_of_plates_to_extract == 4 :

                    # Load Labware
                    tiprack_1 = ctx.load_labware(tipsbox, '8')
                    tiprack_2 = ctx.load_labware(tipsbox, '9')
                    tiprack_3 = ctx.load_labware(tipsbox, '5')
                    tiprack_4 = ctx.load_labware(tipsbox, '6')
                    transfer_tiprack_1 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                    transfer_tiprack_2 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                    transfer_tiprack_3 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                    transfer_tiprack_4 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                    tiprack_9 = ctx.load_labware(tipsbox, protocol_api.OFF_DECK)
                    samples_plate_1 = ctx.load_labware(samples_plate_type, '1')
                    samples_plate_2 = ctx.load_labware(samples_plate_type, '2')
                    samples_plate_3 = ctx.load_labware(samples_plate_type, '3')
                    samples_plate_4 = ctx.load_labware(samples_plate_type, '7')
                    final_plate_1 = ctx.load_labware(final_plate_type, protocol_api.OFF_DECK)
                    final_plate_2 = ctx.load_labware(final_plate_type, protocol_api.OFF_DECK)
                    final_plate_3 = ctx.load_labware(final_plate_type, protocol_api.OFF_DECK)
                    final_plate_4 = ctx.load_labware(final_plate_type, protocol_api.OFF_DECK)
                    reservoir_1 = ctx.load_labware(reservoir_type, '11')
                    water_reservoir = ctx.load_labware(reservoir_type, '10')
                    trash = ctx.fixed_trash['A1']

                    # Load Instrument
                    p300 = ctx.load_instrument('p300_multi_gen2', 'left', tip_racks=[tiprack_1, tiprack_2, tiprack_3, tiprack_4, transfer_tiprack_1, transfer_tiprack_2, transfer_tiprack_3, transfer_tiprack_4, tiprack_9])
                    p300.default_speed = 200

                    # Subseting only the columns with samples in your plates
                    samples_P1 = samples_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                    samples_P2 = samples_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                    samples_P3 = samples_plate_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                    samples_P4 = samples_plate_4.rows()[0][first_column_plate_4 - 1 : last_column_plate_4]
                    plates = [samples_P1, samples_P2, samples_P3, samples_P4]                
                    all_samples = [well for plate in plates for well in plate]

                    mixing_tiprack_1 = tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                    mixing_tiprack_2 = tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                    mixing_tiprack_3 = tiprack_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                    mixing_tiprack_4 = tiprack_4.rows()[0][first_column_plate_4 - 1 : last_column_plate_4]
                    mixing_tipracks = [mixing_tiprack_1, mixing_tiprack_2, mixing_tiprack_3, mixing_tiprack_4]
                    
                    final_plate_01 = final_plate_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                    final_plate_02 = final_plate_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                    final_plate_03 = final_plate_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                    final_plate_04 = final_plate_4.rows()[0][first_column_plate_4 - 1 : last_column_plate_4]
                    final_plates = [final_plate_01, final_plate_02, final_plate_03, final_plate_04]
                    all_final_plates = [well for plate in final_plates for well in plate]

                    transfer_tiprack_01 = transfer_tiprack_1.rows()[0][first_column_plate_1 - 1 : last_column_plate_1]
                    transfer_tiprack_02 = transfer_tiprack_2.rows()[0][first_column_plate_2 - 1 : last_column_plate_2]
                    transfer_tiprack_03 = transfer_tiprack_3.rows()[0][first_column_plate_3 - 1 : last_column_plate_3]
                    transfer_tiprack_04 = transfer_tiprack_4.rows()[0][first_column_plate_4 - 1 : last_column_plate_4]
                    transfer_tipracks = [transfer_tiprack_01, transfer_tiprack_02, transfer_tiprack_03, transfer_tiprack_04]               
                    
                    reservoir_01 = reservoir_1.wells()[0]
                    water_reservoir_01 = water_reservoir.wells()[0]


                    #############################################################################################################
                    #                                           DNA extraction actions
                    #############################################################################################################
                
                    ctx.pause(comment_reagents_1(Metabisulfite, PVPK29, StockLysisSolution_A, StockLysisSolution_B, Sarkosyl, Rnase))
                    ctx.pause(comment_reagents_2(AlcoholIsoamyl, Chloroform, volume_of_water1, volume_of_ethanol95))
                    ctx.pause('''Place full Isopropanol reservoir at -20°C''')
                    ctx.pause(comment_beginning)
                    
                    ctx.pause('''Place Samples plates '1, 2, 3 & 4' on sites 1, 2, 3 & 7''')
                    ctx.pause('''Place tip racks on site 5, 6, 8, 9 and place TE buffer reservoir on site 10''')
                    ctx.pause('''START TE buffer  dispensing''')

                    water_distribute(all_samples, volume = 50, bottom_dispense_heigth = 40, return_tip = True, source =  water_reservoir_01)

                    ctx.pause('''Grind samples on a tyssus-lyser machine''')
                    ctx.pause(comment_spinDown)
                    ctx.pause('''Place the plates back to respective sites 1, 2, 3 & 7. Add Rnase to Extraction buffer and place in reservoir on site 11''')
                    ctx.pause('''START post-grinding Extraction buffer dispensing''')

                    ExtractionBuffer_dispense()
                    
                    ctx.pause('''Seal plates with sealing tape and invert plates 10 times. Spin plates then remove sealing tape and incubate the plates (65C, 60 min)''')
                    ctx.pause('''After incubation, place the plate back to respective sites 1, 2, 3, 7 and place Chloroform reservoir on site 11 and place water reservoir on site 10)''')
                    time_estimation = str(truncate((time_dispensing_chloroform_and_bubbleMixing * total_number_of_columns / 60), 1))
                    ctx.pause(comment_start_Chloro_dispensing(time_estimation))
       
                    if chloroform_buffer_mixing == 'pipette_mixing':
                        dispensing_chloroform_and_pipetteMixing()
                        ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                        ctx.pause('''Remove tip racks on sites 5, 6, 8 & 9 (dispose tips, clean tipracks) and place new tipracks on sites 8 & 9''')
                    if chloroform_buffer_mixing == 'bubble_mixing':
                        dispensing_chloroform_and_bubbleMixing()
                        ctx.pause('''Centrifugate the plate (6000rpm, 10 min).''') 
                        ctx.pause('''Remove tip racks on sites 5, 6, 8 & 9 (dispose tips, clean tipracks) and place new tipracks on sites 8 & 9''')
                    if chloroform_buffer_mixing == 'no_mixing':
                        dispensing_chloroform()
                        ctx.pause('''Mix (vortex carefully) then centrifugate the plate (6000rpm, 10 min).''')
                    
                    ctx.pause('''Remove reservoir on site 10 and place empty plates (96-Deep wells) on sites 4, 5, 6 & 10 (identify the plates and mark site number)''')

                    for tiprack in [tiprack_1, tiprack_2, tiprack_3, tiprack_4]:
                        ctx.move_labware(labware = tiprack, new_location = protocol_api.OFF_DECK)

                    ctx.move_labware(labware = water_reservoir, new_location = protocol_api.OFF_DECK)
                    ctx.move_labware(labware = transfer_tiprack_1, new_location = 8)
                    ctx.move_labware(labware = transfer_tiprack_2, new_location = 9)
                    ctx.move_labware(labware = final_plate_1, new_location = 4)
                    ctx.move_labware(labware = final_plate_2, new_location = 5)
                    ctx.move_labware(labware = final_plate_3, new_location = 6)
                    ctx.move_labware(labware = final_plate_4, new_location = 10)

                    ctx.pause('''When the 10 minutes centrifugation is done place the Sample plates back to respective sites 1, 2, 3 & 7''')

                    # Redefining vectors since we have only two tipracks on the deck
                    plates = [samples_P1, samples_P2]
                    final_plates =[final_plate_01, final_plate_02]
                    transfer_tipracks = [transfer_tiprack_01, transfer_tiprack_02]
                    
                    time_estimation = str(truncate(time_supernatant_transfer * (number_of_columns_plate_1 + number_of_columns_plate_2) / 60, 1))
                    ctx.pause(comment_start_Supernatant_transfer(time_estimation))

                    for samples, final, tiprack in zip(plates, final_plates, transfer_tipracks):
                        Supernatant_transfer(samples, final, tiprack)

                    ctx.pause('''Remove Samples plates on sites 1 & 2 and place new tipracks on sites 1 & 2''')
                
                    for plate in [samples_plate_1, samples_plate_2]:
                        ctx.move_labware(labware = plate, new_location = protocol_api.OFF_DECK)
                
                    ctx.move_labware(labware = transfer_tiprack_3, new_location = 2)
                    ctx.move_labware(labware = transfer_tiprack_4, new_location = 1)

                    time_estimation = str(truncate(time_supernatant_transfer * (number_of_columns_plate_3 + number_of_columns_plate_4) / 60, 1))
                    ctx.pause(comment_start_Supernatant_transfer(time_estimation))

                    # Redefining vectors since we have only two tipracks on the deck
                    plates = [samples_P3, samples_P4]
                    final_plates =[final_plate_03, final_plate_04]
                    transfer_tipracks = [transfer_tiprack_03, transfer_tiprack_04]

                    for samples, final, tiprack in zip(plates, final_plates, transfer_tipracks):
                        Supernatant_transfer(samples, final, tiprack)


                    for plate in [samples_plate_3, samples_plate_4]:
                        ctx.move_labware(labware = plate, new_location = protocol_api.OFF_DECK)
                    
                    # Redefining plate vectors
                    final_plates =[final_plate_01, final_plate_02, final_plate_03, final_plate_04]
                    transfer_tipracks = [transfer_tiprack_01, transfer_tiprack_02, transfer_tiprack_03, transfer_tiprack_04]

                    ctx.pause('''Remove plates on sites 3 & 7 and Place a new tiprack on site 3''')
                    ctx.move_labware(labware = tiprack_9, new_location = 3)

                    ctx.pause('''Remove Chloroform reservoir on site 11 and place cold Isopropanol reservoir on site 11''')
                    ctx.pause('''START Isopropanol dispensing to plates on sites 4, 5, 6 & 10''')

                    isopropanol_dispensing()

                    ctx.pause('''At this point the plates (sites 4, 5, 6 & 10) can be sealed and stored at -20C and the DNA extraction completed the next day (or later)'''  )
                    ctx.pause('''Seal and invert 20x the plates on sites 4, 5, 6 & 10 then centrifugate (10min, 6000 rpm)''') # and prepare cold 70 percent ethanol

                    # Evaluating if the Isopropanol discarding is done by plate inversion or py pipetting off
                    if pipetteOff_isopropanol == True:
                        ctx.pause('''When centrifugation is done, place the plate back to respective sites 4, 5, 6 & 10''')
                        time_estimation = str(truncate(time_isopropanol_discarding * total_number_of_columns / 60, 1))
                        ctx.pause(comment_start_Isopropanol_discarding(time_estimation))

                        for final, tiprack in zip(final_plates, transfer_tipracks):
                            isopropanol_discarding(final, tiprack)

                    if pipetteOff_isopropanol == False:
                        ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates 10s''')
                        ctx.pause('''Place plates back to respective sites 4, 5, 6 & 10''')

                    ctx.pause('''Remove Isopropanol reservoir on site 11 and place 70'%'Ethanol reservoir to site 11''')
                    ctx.pause('''START Ethanol dispensing to plates on site 4, 5, 6 & 10''')

                    ethanol_dispensing()

                    ctx.pause('''Centrifugate the plates for 10 minutes at 6000 rpm''')

                    # Ethanol discarding method (by inversion or py pipetting off)       
                    if pipetteOff_ethanol == True:
                        ctx.pause('''When centrifugation is done, place the plates back to respective sites 4, 5, 6 & 10''')
                        time_estimation = str(truncate(time_ethanol_discarding * total_number_of_columns / 60, 1))
                        ctx.pause(comment_start_Ethanol_discarding(time_estimation))
                        for final, tiprack in zip(final_plates, transfer_tipracks):
                            ethanol_discarding(final, tiprack)
                    if pipetteOff_ethanol == False:
                        ctx.pause('''Gently invert the plates to poor off the supernatant then centrifuge the plates (6000rpm, 10sec)''')

                    ctx.pause('''Evaporate ethanol (20min, 45C). Remove and clean tipracks on sites 1, 2, 8 & 9 (place tips in trash)''')
                    ctx.pause('''Prepare for Elution buffer dispensing: remove Ethanol reservoir on site 11 and place Elution buffer on site 11''')
                    ctx.pause('''When evaporation is done, place plates back to respective sites 4, 5, 6 and 10''')
                    ctx.pause('''START Elution buffer dispensing''')

                    ctx.comment('\n~~~~~~~~~~~~~~ Dispensing Elution buffer to samples ~~~~~~~~~~~~~~')

                    EB_distribute(all_final_plates, volume = elution_buffer_volume, bottom_dispense_heigth = 16, return_tip = False, source = reservoir_01)
                    ctx.pause('''Centrifuge plates at 4000rpm 5s. Seal plates with tape and store plates at 4C. DNA extraction completed''')
                    ctx.comment('\n~~~~~~~~~~~~~~Protocol Complete~~~~~~~~~~~~~~\n')
