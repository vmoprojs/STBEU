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
    //checkError(err, "getting platforms");
    
    // Enumerate devices
    unsigned numDevices = 0;
    for (int i = 0; i < numPlatforms; i++)
    {
        cl_uint num = 0;
        err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                             MAX_DEVICES-numDevices, devices+numDevices, &num);
        //checkError(err, "getting deviceS");
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
                Rprintf("No devices found.\n");
            }
            else
            {
                Rprintf("\n");
                Rprintf("Devices:\n");
                for (int i = 0; i < numDevices; i++)
                {
                    char name[MAX_INFO_STRING];
                    getDeviceName(devices[i], name);
                    Rprintf("%2d: %s\n", i, name);
                }
                Rprintf("\n");
            }
            //exit(0);
        }
        else if (!strcmp(argv[i], "--device"))
        {
            if (++i >= argc || !parseUInt(argv[i], deviceIndex))
            {
                Rprintf("Invalid device index\n");
                //exit(1);
            }
        }
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
        {
            Rprintf("\n");
            Rprintf("Usage: ./program [OPTIONS]\n\n");
            Rprintf("Options:\n");
            Rprintf("  -h  --help               Print the message\n");
            Rprintf("      --list               List available devices\n");
            Rprintf("      --device     INDEX   Select device at INDEX\n");
            Rprintf("\n");
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
 *-----------------------------------------------------------------------*/


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
        Rprintf("Error during operation '%s', ", operation);
        Rprintf("in '%s' on line %d\n", filename, line);
        Rprintf("Error code was \"%s\" (%d)\n", err_code(err), err);
    }
}


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
    char vendor[1024];                //this strirng will hold a platforms vendor
    
    
    for (i = 0; i < platformCount; i++) {
        
        Rprintf("Platform:\t\t%u\n\n", i);
        clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(vendor), vendor, NULL);
        Rprintf("\tPlatform Vendor:\t%s\n", vendor);
        
        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
        
        //int gpu[deviceCount];
        //int countergpu=0;
        SEP;
        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
            SEP;
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
            SEP;
            
        }
        SEP;
        free(devices);
        
    }
    
    free(platforms);
    return 0;
    
}



// ================================ Start  Create Binary Kernel

void create_binary_kernel(int *dev, char **fname)
{
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
    char CL[5];
    char f_nameCL[100];
    strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    //strcat(f_nameCL, CL);
    //Rprintf("ARCHIVO CL  %s\n", f_nameCL);
    
    
    char *kernelsource = getKernelSource(f_nameCL);
    
    program = clCreateProgramWithSource(context,1,(const char **) &kernelsource, NULL, &err);
    if(err!=CL_SUCCESS){Rprintf("Failed clCreateProgramWithSource\n");}
    
    err = clCompileProgram(program, 1, &device, "-I ./", 0, NULL, NULL, NULL, NULL);
    //err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*900];
        
        Rprintf("CREATE Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("CREATE Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("CREATE Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("CREATE Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }

    FILE *f;
    char *binary;
    size_t binary_size;
    
    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binary_size, NULL);
    
    binary = malloc((binary_size)); // FUCK!!!
    clGetProgramInfo(program, CL_PROGRAM_BINARIES, binary_size, &binary, NULL);
    //if(err!=CL_SUCCESS){RRprintf("Failed to get CL_PROGRAM_BINARIES\n");}
    
    strcat(*fname, ".bin");
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


void create_binary_kernel_GPU(int *dev, char **fname)
{
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
    char CL[5];
    char f_nameCL[100];
    strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    //strcat(f_nameCL, CL);
    //Rprintf("ARCHIVO CL  %s\n", f_nameCL);
    
    
    char *kernelsource = getKernelSource(f_nameCL);
    
    program = clCreateProgramWithSource(context,1,(const char **) &kernelsource, NULL, &err);
    if(err!=CL_SUCCESS){Rprintf("Failed clCreateProgramWithSource\n");}
    
    //err = clCompileProgram(program, 1, &device, "-I ./", 0, NULL, NULL, NULL, NULL);
    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*900];
        
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

    FILE *f;
    char *binary;
    size_t binary_size;
    
    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binary_size, NULL);
    
    binary = malloc((binary_size)); // FUCK!!!
    clGetProgramInfo(program, CL_PROGRAM_BINARIES, binary_size, &binary, NULL);
    //if(err!=CL_SUCCESS){RRprintf("Failed to get CL_PROGRAM_BINARIES\n");}
    
    strcat(*fname, ".bin");
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
        //fRprintf(stderr, "Error: Could not open kernel source file\n");
        Rprintf("Error: Could not open kernel source file\n");
        //exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);
    
    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        //fRprintf(stderr, "Error: Could not allocate memory for source string\n");
        //exit(EXIT_FAILURE);
        Rprintf("Error: Could not allocate memory for source string\n");
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

/**************** sum total *********/
double sum_total(double *arr, int ngrid)
{
    double sol = 0.0;
    for (int i=0; i<ngrid; i++)
    {
        sol += arr[i];
    }
    return sol;
}

void param_OCL(int *npts,int *ntime,int *flagcor,int *flagnuis,int *npar,int *weigthed, int *dist,double *maxtime,double *maxdist,double *parcor,double *nuis,int *nparc,int *cormod,int *int_par, double *dou_par,int *nparnuis)
{
    int_par[0] = npts[0];
    int_par[1] = ntime[0];
    int_par[2] = flagcor[0];
    int_par[3] = flagcor[1];
    int_par[4] = flagcor[2];
    int_par[5] = flagcor[3];
    int_par[6] = flagcor[4];
    int_par[7] = flagnuis[0];
    int_par[8] = flagnuis[1];
    int_par[9] = flagnuis[2];
    int_par[10] = npar[0];//npar
    int_par[11] = weigthed[0];
    int_par[12] = dist[0];
    int_par[13] = nparc[0];
    int_par[14] = nparnuis[0];
    int_par[15] = cormod[0];
    
    int_par[16] = flagcor[5];
    int_par[17] = flagcor[6];
    
    dou_par[0] = maxtime[0];
    dou_par[1] = maxdist[0];
    dou_par[2] = parcor[0];
    dou_par[3] = parcor[1];//parcor
    dou_par[4] = parcor[2];//parcor
    dou_par[5] = parcor[3];//parcor
    dou_par[6] = parcor[4];//parcor
    dou_par[7] = nuis[0];
    dou_par[8] = nuis[1];
    dou_par[9] = nuis[2];
    
    dou_par[10] = parcor[5];//parcor
    dou_par[11] = parcor[6];//parcor
    
    
}


void exec_kernel(int *npts, int *ntime,double *h_t, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *h_data,int *weigthed, double *mom_cond, int *dist, double *h_x,double *h_y,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev,char **fname, int *nparc, int *nparnuis)
{
    
   int i;
    
    int *int_par;
        double *dou_par;
        int_par = (int*)Calloc((64), int *);
        dou_par = (double*)Calloc((64), double *);
    
    
    param_OCL(npts,ntime,flagcor,flagnuis,npar,weigthed,dist,maxtime,maxdist,parcor,nuis,nparc,cormod,int_par,dou_par,nparnuis);
    
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    
    char name[MAX_INFO_STRING];
    getDeviceName(device, name);
    
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
   
    
    FILE *fp;
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    char f_nameCL[100];
    //strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    
     //strcat(*fname, ".bin");
    strcat(f_nameCL, ".bin");
    fp = fopen(f_nameCL, "r");//"DouExp.cl.bin"
    
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    free(binary_buf);
    checkError(err, "Creating Creating Binary Program");
    // Build the program
    err = clBuildProgram(program,1, &device, NULL, NULL, &err);
     //checkError(err, "Creating Building Program");
    if (err != CL_SUCCESS)
       {
           size_t len;
           char buffer[2048*900];
           
           Rprintf("EXEC Error: Failed to build program executable!\n%s\n", err_code(err));
           clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
           SEP;
           Rprintf("EXEC Build Log:\n%s\n", buffer);
           SEP;
           clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
           Rprintf("EXEC Build Status:\n%s\n", buffer);
           SEP;
           clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
           Rprintf("EXEC Build Options:\n%s\n", buffer);
           SEP;
           //return EXIT_FAILURE;
       }
    
    
    
    // Create the compute kernel from the program
    kernel = clCreateKernel(program, "scalarspaceocl", &err);
    checkError(err, "Creating kernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");
   
    
   
    // Creating buffers
    size_t coords_buff = sizeof(double) * int_par[0];
    size_t coordt_buff = sizeof(double) * int_par[1];
    size_t data_buff = sizeof(double) * (int_par[0]*int_par[1]);
    size_t int_par_buff = sizeof(int) * ((int)64);
    size_t dou_par_buff = sizeof(double) *((double)64);
    
    size_t g1,g2;
    const int ll1 =local_wi[0];
    const int ll2 =local_wi[1];
    g1 = int_par[0] + (ll1 - (int_par[0] & (ll1-1))); // SPACE
    g2 = int_par[1] + (ll2 - (int_par[1] & (ll2-1))); //TIME
    
    
    //printf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local[2] = {ll1,ll2};
    size_t global[2] = {g1,g2};
    
    int length = g1*g2;
    size_t length_buff = sizeof(double)* (length);
        
    double *h_mom_cond0,*h_mom_cond1,*h_mom_cond2,*h_mom_cond3;
    h_mom_cond0= (double*)calloc(length, sizeof(double));
    h_mom_cond1= (double*)calloc(length, sizeof(double));
    h_mom_cond2= (double*)calloc(length, sizeof(double));
    h_mom_cond3= (double*)calloc(length, sizeof(double));
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    checkError(err, "Creating buffer device X");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, coords_buff, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device X");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, coords_buff, (void*)h_y, 0, NULL, NULL);
    //checkError(err, "Creating buffer device Y");
    
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, data_buff, (void*)h_data, 0, NULL, NULL);
    //checkError(err, "Creating buffer device data");
    
    cl_mem d_t = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_t, CL_TRUE, 0, coordt_buff, (void*)h_t, 0, NULL, NULL);
    //checkError(err, "Creating buffer device time");
    
    cl_mem m0_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    //checkError(err, "Creating buffer device m0_sol");
    err = clEnqueueWriteBuffer(commands, m0_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond0, 0, NULL, NULL);
    //checkError(err, "Writing buffer device m0_sol");
    
    
    cl_mem m1_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    //checkError(err, "Creating buffer device m1_sol");
    err = clEnqueueWriteBuffer(commands, m1_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond1, 0, NULL, NULL);
    //checkError(err, "Writing buffer device m1_sol");
    
    cl_mem m2_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    //checkError(err, "Creating buffer device m2_sol");
    err = clEnqueueWriteBuffer(commands, m2_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond2, 0, NULL, NULL);
   // checkError(err, "Writing buffer device m2_sol");
    
    cl_mem m3_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    //checkError(err, "Creating buffer device m3_sol");
    err = clEnqueueWriteBuffer(commands, m3_sol, CL_TRUE, 0, length_buff, (void*)h_mom_cond3, 0, NULL, NULL);
    //checkError(err, "Writing buffer device m3_sol");
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    //checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    //checkError(err, "Creating buffer device double params");
    //Rprintf("??????INSIDE FUNCTION\n");
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_t); //coordt
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_x); //coordx
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_y); //coordx
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data); //data
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &m0_sol); //mom_cond0
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &m1_sol); //mom_cond1
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &m2_sol); //mom_cond2
    err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &m3_sol); //mom_cond3
    err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_int_par); //mom_cond3
    err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_dou_par); //mom_cond3
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global,local, 0, NULL, NULL);
    clFinish(commands);
    err = clEnqueueReadBuffer(commands, m0_sol, CL_TRUE, 0,length_buff, h_mom_cond0, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m1_sol, CL_TRUE, 0,length_buff, h_mom_cond1, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m2_sol, CL_TRUE, 0,length_buff, h_mom_cond2, 0, NULL, NULL);
    err = clEnqueueReadBuffer(commands, m3_sol, CL_TRUE, 0,length_buff, h_mom_cond3, 0, NULL, NULL);
    clFinish(commands);
    
    double m0=0,m1=0,m2=0,m3=0;
    for(i=0; i<(int_par[0]*int_par[1]); i++)
    {
        m0 += h_mom_cond0[i];
        m1 += h_mom_cond1[i];
        m2 += h_mom_cond2[i];
        m3 += h_mom_cond3[i];
        /*if(isnan(m0) || isnan(m1) || isnan(m2) || isnan(m3))
        {
            //printf("npts[0] %d, ntime[0]  %d,    flagcor[0]  %d,    flagcor[1]  %d,   flagcor[2]  %d,    flagcor[3]  %d,   flagcor[4]  %d,    flagnuis[0]  %d,   flagnuis[1]  %d,   flagnuis[2]  %d,    npar[0]  %d,  weigthed[0]  %d,   dist[0]  %d,   nparc[0]  %d,   nparnuis[0]  %d,    cormod[0]  %d\n\n",npts[0], ntime[0],    flagcor[0],    flagcor[1],   flagcor[2],    flagcor[3],   flagcor[4],    flagnuis[0],   flagnuis[1],   flagnuis[2],    npar[0] ,  weigthed[0],   dist[0],   nparc[0],   nparnuis[0],    cormod[0]);
            
           //printf("maxtime[0] %f,    maxdist[0] %f,    parcor[0] %f,    parcor[1] %f,    parcor[2] %f,    parcor[3] %f,    parcor[4] %f,    nuis[0] %f,    nuis[1] %f,    nuis[2] %f,        parcor[5] %f,    parcor[6] %f\n\n",maxtime[0],    maxdist[0],    parcor[0],    parcor[1],    parcor[2],    parcor[3],    parcor[4],    nuis[0],    nuis[1],    nuis[2],        parcor[5],    parcor[6]);
            
            
            printf("\nm0: %f m1: %f m2: %f m3: %f \n",m0,m1,m2,m3);
            printf("\nHm0: %f Hm1: %f Hm2: %f Hm3: %f \n",h_mom_cond0[i],h_mom_cond1[i],h_mom_cond2[i],h_mom_cond3[i]);
        }*/
        
    }
    mom_cond[0] =m0;mom_cond[1] =m1;mom_cond[2] =m2;mom_cond[3] =m3;
    // clean up inside kernels

    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    //clReleaseDevice(device);
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_t);
    clReleaseMemObject(m0_sol);
    clReleaseMemObject(m1_sol);
    clReleaseMemObject(m2_sol);
    clReleaseMemObject(m3_sol);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
    //clReleaseMemObject(err);
    free(h_mom_cond0);
    free(h_mom_cond1);
    free(h_mom_cond2);
    free(h_mom_cond3);
}



