################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../CenData.o \
../Density.o \
../MyMath.o \
../OPT.o \
../PHI.o \
../WKM.o \
../main.o 

CPP_SRCS += \
../CenData.cpp \
../Density.cpp \
../MyMath.cpp \
../OPT.cpp \
../PHI.cpp \
../WKM.cpp \
../main.cpp 

OBJS += \
./CenData.o \
./Density.o \
./MyMath.o \
./OPT.o \
./PHI.o \
./WKM.o \
./main.o 

CPP_DEPS += \
./CenData.d \
./Density.d \
./MyMath.d \
./OPT.d \
./PHI.d \
./WKM.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


