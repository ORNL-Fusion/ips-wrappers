	#-------------------------------------------------------------------------------
#
#  This utility writes messages to the console based on a predfined priority 
#  level. The uses can set the amount of screen output in the work flow config
#  file. Levels are
#
#  quiet   - Write output reguardless of the priority level.
#  verbose - Only write output if the config file output level is set to 
#            verbose.
#
#  Users can choose the level of screen output by setting the OUTPUT_LEVEL 
#  variable in the config file. If the variable is not set, assume verbose 
#  output.
#
#  A developer can set the priority of output messages using the same levels.
#
#-------------------------------------------------------------------------------
def screen_output(component, level, message):
    try:
        output_level = component.services.get_config_param('OUTPUT_LEVEL')

        if output_level == 'verbose':
            print(message)
        elif level == 'quiet':
            print(message)
        elif level != 'verbose':
            raise Exception

    except KeyError:
        if level == 'verbose' or level == 'quiet':
            print(message)
        else:
            raise Exception
