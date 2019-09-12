#include "header.h"

/*------------------------------------------------------------------------------
 *
 * Name:       device_picker.h
 *
 * Purpose:    Provide a simple CLI to specify an OpenCL device at runtime
 *
 * Note:       Must be included AFTER the relevant OpenCL header
 *             See one of the Matrix Multiply exercises for usage
 *
 * HISTORY:    Method written by James Price, October 2014
 *             Extracted to a common header by Tom Deakin, November 2014
 */




unsigned getDeviceList(cl_device_id devices[MAX_DEVICES])
{
    cl_int err;
    
    // Get list of platforms
    cl_uint numPlatforms = 0;
    cl_platform_id platforms[MAX_PLATFORMS];
    err = clGetPlatformIDs(MAX_PLATFORMS, platforms, &numPlatforms);
    checkError(err, "getting platforms");
    
    // Enumerate devices
    unsigned numDevices = 0;
    for (int i = 0; i < numPlatforms; i++)
    {
        cl_uint num = 0;
        err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                             MAX_DEVICES-numDevices, devices+numDevices, &num);
        checkError(err, "getting deviceS");
        numDevices += num;
    }
    
    return numDevices;
}

void getDeviceName(cl_device_id device, char name[MAX_INFO_STRING])
{
    cl_device_info info = CL_DEVICE_NAME;
    
    // Special case for AMD
#ifdef CL_DEVICE_BOARD_NAME_AMD
    clGetDeviceInfo(device, CL_DEVICE_VENDOR, MAX_INFO_STRING, name, NULL);
    if (strstr(name, "Advanced Micro Devices"))
        info = CL_DEVICE_BOARD_NAME_AMD;
#endif
    
    clGetDeviceInfo(device, info, MAX_INFO_STRING, name, NULL);
}


int parseUInt(const char *str, cl_uint *output)
{
    char *next;
    *output = strtoul(str, &next, 10);
    return !strlen(next);
}

void parseArguments(int argc, char *argv[], cl_uint *deviceIndex)
{
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--list"))
        {
            // Get list of devices
            cl_device_id devices[MAX_DEVICES];
            unsigned numDevices = getDeviceList(devices);
            
            // Print device names
            if (numDevices == 0)
            {
                printf("No devices found.\n");
            }
            else
            {
                printf("\n");
                printf("Devices:\n");
                for (int i = 0; i < numDevices; i++)
                {
                    char name[MAX_INFO_STRING];
                    getDeviceName(devices[i], name);
                    printf("%2d: %s\n", i, name);
                }
                printf("\n");
            }
            //exit(0);
        }
        else if (!strcmp(argv[i], "--device"))
        {
            if (++i >= argc || !parseUInt(argv[i], deviceIndex))
            {
                printf("Invalid device index\n");
                //exit(1);
            }
        }
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
        {
            printf("\n");
            printf("Usage: ./program [OPTIONS]\n\n");
            printf("Options:\n");
            printf("  -h  --help               Print the message\n");
            printf("      --list               List available devices\n");
            printf("      --device     INDEX   Select device at INDEX\n");
            printf("\n");
            //exit(0);
        }
    }
}










//#pragma once
/*----------------------------------------------------------------------------
 *
 * Name:     err_code()
 *
 * Purpose:  Function to output descriptions of errors for an input error code
 *           and quit a program on an error with a user message
 *
 *
 * RETURN:   echoes the input error code / echos user message and exits
 *
 * HISTORY:  Written by Tim Mattson, June 2010
 *           This version automatically produced by genErrCode.py
 *           script written by Tom Deakin, August 2013
 *           Modified by Bruce Merry, March 2014
 *           Updated by Tom Deakin, October 2014
 *               Included the checkError function written by
 *               James Price and Simon McIntosh-Smith
 *
 *----------------------------------------------------------------------------
 */
/*
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#ifdef __cplusplus
#include <cstdio>
#endif
 */

const char *err_code (cl_int err_in)
{
    switch (err_in) {
        case CL_SUCCESS:
            return (char*)"CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:
            return (char*)"CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE:
            return (char*)"CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:
            return (char*)"CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            return (char*)"CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:
            return (char*)"CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:
            return (char*)"CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            return (char*)"CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
            return (char*)"CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
            return (char*)"CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            return (char*)"CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
            return (char*)"CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE:
            return (char*)"CL_MAP_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            return (char*)"CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            return (char*)"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case CL_INVALID_VALUE:
            return (char*)"CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:
            return (char*)"CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:
            return (char*)"CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:
            return (char*)"CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:
            return (char*)"CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:
            return (char*)"CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:
            return (char*)"CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:
            return (char*)"CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:
            return (char*)"CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            return (char*)"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:
            return (char*)"CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:
            return (char*)"CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:
            return (char*)"CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:
            return (char*)"CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:
            return (char*)"CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:
            return (char*)"CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:
            return (char*)"CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:
            return (char*)"CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
            return (char*)"CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:
            return (char*)"CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
            return (char*)"CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
            return (char*)"CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
            return (char*)"CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:
            return (char*)"CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
            return (char*)"CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:
            return (char*)"CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET:
            return (char*)"CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:
            return (char*)"CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:
            return (char*)"CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:
            return (char*)"CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:
            return (char*)"CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:
            return (char*)"CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:
            return (char*)"CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE:
            return (char*)"CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_PROPERTY:
            return (char*)"CL_INVALID_PROPERTY";
            
        default:
            return (char*)"UNKNOWN ERROR";
    }
}


void check_error(cl_int err, const char *operation, char *filename, int line)
{
    if (err != CL_SUCCESS)
    {
        //fprintf(stderr, "Error during operation '%s', ", operation);
        //fprintf(stderr, "in '%s' on line %d\n", filename, line);
        //fprintf(stderr, "Error code was \"%s\" (%d)\n", err_code(err), err);
        //exit(EXIT_FAILURE);
        printf("Error during operation '%s', \n", operation);
    }
}


#define checkError(E, S) check_error(E,S,__FILE__,__LINE__)


// ************* **************       ***************   OPENCL
int DevOpenCL()
{
    
    int i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    cl_device_fp_config fp;
    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount);
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
    clGetPlatformIDs(platformCount, platforms, NULL);
    cl_device_type dt;
    
    
    for (i = 0; i < platformCount; i++) {
        
        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
        
        //int gpu[deviceCount];
        //int countergpu=0;
        Rprintf("-------------------------------------------------------------\n");
        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
            Rprintf("-------------------------------------------------------------\n");
            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            Rprintf("%d.\tCL_DEVICE_NAME\tDevice: \t%s\n", j, value);
            free(value);
            
            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_VERSION\tHardware version: \t%s\n", j, 1, value);
            free(value);
            
            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DRIVER_VERSION\tSoftware version: \t%s\n", j, 2, value);
            free(value);
            
            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_OPENCL_C_VERSION\tOpenCL C version: \t%s\n", j, 3, value);
            free(value);
            
            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_COMPUTE_UNITS\tParallel compute units: \t%d\n", j, 4, maxComputeUnits);
            
            
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_SIZE\tGlobal Memory (MB):\t%llu\n",j, 5,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_CACHE_SIZE\tGlobal Memory Cache (MB):\t%llu\n",j, 6,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_LOCAL_MEM_SIZE\tLocal Memory (KB):\t%llu\n",j, 7,long_entries/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_CLOCK_FREQUENCY\tMax clock (MHz) :\t%llu\n",j, 8,long_entries);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_GROUP_SIZE\tMax Work Group Size:\t%zu\n",j, 9,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_SIZES\tMax Work Item Size:\t%zu\n",j,10,p_size);
            clGetDeviceInfo(devices[j],CL_KERNEL_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_KERNEL_WORK_GROUP_SIZE\tMax kernel Work group Size:\t%zu\n",j, 11,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_DIMENSIONS\tMax Dev dim:\t%zu\n",j, 12,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_MEM_ALLOC_SIZE\tMax Buffer size (Mb):\t%zu\n",j, 13,p_size/1000000);
            clGetDeviceInfo(devices[j],CL_DEVICE_DOUBLE_FP_CONFIG,sizeof(cl_device_fp_config),&fp,NULL);
            Rprintf("%d.%d\tSupports double precision floating-point? %s\n",j, 14,fp != 0 ? "yes" : "no");
            clGetDeviceInfo(devices[j],CL_DEVICE_TYPE,sizeof(cl_device_type),&dt,NULL);
            Rprintf("%d.%d\tCL_DEVICE_TYPE: %s\n",j, 15,dt & CL_DEVICE_TYPE_GPU ? "GPU" : "CPU");
            //clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CONSTANT_ARGS,sizeof(size_t),&p_size,NULL);
            //Rprintf("%d.%d\tCL_DEVICE_MAX_CONSTANT_ARGS\tMax kernel args:\t%zu\n",j, 15,p_size);
            Rprintf("-------------------------------------------------------------\n");
            
        }
        Rprintf("-------------------------------------------------------------\n");
        free(devices);
        
    }
    
    free(platforms);
    return 0;
    
}





// ================================ Start  Create Binary Kernel

void create_binary_kernel(int *dev, char **fname)
{
    //printf("ARCHIVO fname  %s\n", *fname);
    // Context, program, build:
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    //cl_kernel kernel;
    
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list') Compilation!\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    // Create a compute context
    /*cl_platform_id platform=NULL;
     cl_context_properties ctx_properties[] = {
     CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0
     };*/
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    
    
    //fclose(fp);*/
    //char CL[5];
    char f_nameCL[100];
    //strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    //strcat(f_nameCL, CL);
    //printf("ARCHIVO CL  %s\n", f_nameCL);
    
    
    char *kernelsource = getKernelSource(f_nameCL);
    
    program = clCreateProgramWithSource(context,1,(const char **) &kernelsource, NULL, &err);
    if(err!=CL_SUCCESS){Rprintf("Failed clCreateProgramWithSource\n");}
    
    err = clCompileProgram(program, 1, &device, "-I ./", 0, NULL, NULL, NULL, NULL);
    //err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*200];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    //kernel = clCreateKernel(program, *fname, &err);
    
    
    //size_t work_group_size;
    //err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    //checkError(err, "Getting kernel work group info");
    //Rprintf("Recommended Local: %lu\n", work_group_size);
    
    
    FILE *f;
    char *binary;
    size_t binary_size;
    
    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binary_size, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARY_SIZES\n");}
    binary = malloc(binary_size);
    clGetProgramInfo(program, CL_PROGRAM_BINARIES, binary_size, &binary, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARIES\n");}
    
    strcat(*fname, ".bin");
    // printf("ARCHIVO CL  B %s\n", *fname);
    f = fopen(*fname, "w");
    fwrite(binary, binary_size, 1, f);
    fclose(f);
    
    
    /* Finalization */
    clFlush(commands);
    clFinish(commands);
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    //clReleaseKernel(kernel);
    Rprintf("Successful Binary Compilation!!!\n");
}

// ================================ End  Create Binary Kernel



// OPEN CL ///

char * getKernelSource(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        //fprintf(stderr, "Error: Could not open kernel source file\n");
        printf("Error: Could not open kernel source file\n");
        //exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);
    
    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        //fprintf(stderr, "Error: Could not allocate memory for source string\n");
        //exit(EXIT_FAILURE);
        printf("Error: Could not allocate memory for source string\n");
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

/**************** sum total *********/
float sum_total(float *arr, int ngrid)
{
    float sol = 0.0;
    for (int i=0; i<ngrid; i++)
    {
        sol += arr[i];
    }
    return sol;
}
