import string

cql3d_output_file_tmp = "'ITER_sn2                                                                                                                                                                                                                                       '"
cql3d_output_file_tmp = cql3d_output_file_tmp.replace(',','').replace("'","")
cql3d_output_file_tmp = cql3d_output_file_tmp.rstrip()
cql3d_output_file = cql3d_output_file_tmp + ".nc"
       
print(cql3d_output_file)       
