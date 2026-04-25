#######################################
########### Dependencies ##############
#######################################

$(OBJ_DIR)/postprocessing/src/energy.o: $(OBJ_DIR)/physical/src/hamiltonians.o \
																				 $(OBJ_DIR)/physical/src/parameters.o \
																				 $(OBJ_DIR)/physical/src/utilities.o \
																				 $(OBJ_DIR)/input_output/src/writers.o \
																				 $(OBJ_DIR)/input_output/src/reader.o \
																				 $(OBJ_DIR)/integrate/src/local_integrand.o \
																				 $(OBJ_DIR)/self_consistency/self_consistency.o \
																				 $(OBJ_DIR)/input_output/src/logger.o \
																				 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/postprocessing/src/symmetry.o: $(OBJ_DIR)/physical/src/hamiltonians.o \
																					 $(OBJ_DIR)/physical/src/parameters.o \
																					 $(OBJ_DIR)/physical/src/utilities.o \
																					 $(OBJ_DIR)/input_output/src/writers.o \
																					 $(OBJ_DIR)/input_output/src/reader.o \
																					 $(OBJ_DIR)/integrate/src/local_integrand.o \
																					 $(OBJ_DIR)/self_consistency/self_consistency.o \
																					 $(OBJ_DIR)/input_output/src/logger.o \
																					 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/postprocessing/src/topology.o: $(OBJ_DIR)/physical/src/hamiltonians.o \
																					 $(OBJ_DIR)/physical/src/parameters.o \
																					 $(OBJ_DIR)/physical/src/utilities.o \
																					 $(OBJ_DIR)/input_output/src/writers.o \
																					 $(OBJ_DIR)/input_output/src/reader.o \
																					 $(OBJ_DIR)/integrate/src/local_integrand.o \
																					 $(OBJ_DIR)/self_consistency/self_consistency.o \
																					 $(OBJ_DIR)/input_output/src/logger.o \
																					 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/physical/src/utilities.o: $(OBJ_DIR)/physical/src/parameters.o \
																			$(OBJ_DIR)/types/types.o

$(OBJ_DIR)/physical/src/parameters.o:

$(OBJ_DIR)/types/types.o: $(OBJ_DIR)/physical/src/parameters.o

$(OBJ_DIR)/physical/src/hamiltonians.o: $(OBJ_DIR)/physical/src/utilities.o \
																				 $(OBJ_DIR)/physical/src/parameters.o \
																				 $(OBJ_DIR)/input_output/src/reader.o \
																				 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/input_output/src/reader.o: $(OBJ_DIR)/physical/src/parameters.o \
																			 $(OBJ_DIR)/input_output/src/logger.o \
																			 $(OBJ_DIR)/types/types.o \
																			 $(OBJ_DIR)/physical/src/utilities.o \
																			 $(OBJ_DIR)/physical/src/interaction_factory.o

$(OBJ_DIR)/input_output/src/logger.o:

$(OBJ_DIR)/input_output/src/writers.o: $(OBJ_DIR)/physical/src/parameters.o \
																				$(OBJ_DIR)/input_output/src/reader.o \
																				$(OBJ_DIR)/types/types.o

$(OBJ_DIR)/integrate/src/local_integrand.o: $(OBJ_DIR)/physical/src/parameters.o \
																						 $(OBJ_DIR)/physical/src/utilities.o \
																						 $(OBJ_DIR)/physical/src/hamiltonians.o \
																						 $(OBJ_DIR)/input_output/src/writers.o \
																						 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/self_consistency/self_consistency.o: $(OBJ_DIR)/physical/src/parameters.o \
																								 $(OBJ_DIR)/input_output/src/reader.o \
																								 $(OBJ_DIR)/input_output/src/logger.o \
																								 $(OBJ_DIR)/types/types.o \
																								 $(OBJ_DIR)/input_output/src/writers.o

$(OBJ_DIR)/physical/src/interaction_factory.o: $(OBJ_DIR)/physical/src/parameters.o \
																								$(OBJ_DIR)/physical/src/utilities.o \
																								$(OBJ_DIR)/types/types.o \
																								$(OBJ_DIR)/input_output/src/logger.o
