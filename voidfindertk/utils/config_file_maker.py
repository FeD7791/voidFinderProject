from configparser import ConfigParser

def config_file_maker(element_name,file_buffer,**kwargs):
    '''
    path: place where the input file is (always has vars.conf as name)
    '''
    config_element = ConfigParser()
    config_element.optionxform = str
    config_element[f"{element_name}"] = kwargs
    config_element.write(file_buffer)
