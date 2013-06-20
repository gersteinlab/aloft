// To compile on Louise:
// gcc gerprate.c -IPython-2.7.3/ -IPython-2.7.3/Include/ -LPython-2.7.3/ -l python2.7 -lpthread -lm -ldl -lutil -fPIC -shared -o gerprate.so
// For platforms where python-config exists, you can get the list of libraries to link via python-config --libs

// when building python, you want to pass CFLAGS=-fPIC to configure

#ifdef __APPLE__
    #include <Python/Python.h>
#else
    #include "Python.h"
#endif

static float *list = NULL;
static FILE *cachedFile = NULL;
static char isByteOrderLittleEndian = 0;

// Returns 1 if the machine is little endian, otherwise 0 if big endian
static char isLittleEndian(void)
{
    int testInteger = 1;
    return (*((unsigned char *)(&testInteger)) == 1);
}

static PyObject *freeMemory(PyObject *self)
{
    if (list)
    {
        free(list);
        list = NULL;
    }
    
    if (cachedFile)
    {
        fclose(cachedFile);
        cachedFile = NULL;
    }
    
    return Py_BuildValue("");
}

static PyObject *buildList(PyObject *self, PyObject *args)
{
    freeMemory(self);
    
    const char *filepath = NULL;
    const char *outfilepath = NULL;
    if (!PyArg_ParseTuple(args, "ss", &filepath, &outfilepath))
    {
        printf("Failed to parse arguments\n");
        return NULL;
    }
    
    cachedFile = fopen(outfilepath, "r");
    if (cachedFile)
    {
        printf("Opened cache for %s\n", outfilepath);
        return Py_BuildValue("");
    }
    
    FILE *file = fopen(filepath, "r");
    if (!file)
    {
        printf("Failed to open file %s\n", filepath);
        return NULL;
    }
    
    unsigned int listMaxCount = 1000000000;
    
    list = malloc(sizeof(float) * listMaxCount);
    if (!list)
    {
        printf("Failed to allocate enough memory for list\n");
        fclose(file);
        return NULL;
    }
    
    // start at index 1 to match line numbers
    list[0] = 0.0;
    unsigned int listIndex = 1;
    
    // read file in chunks
    char buffer[65536];
    size_t numberOfBytesRead = 0;
    while ((numberOfBytesRead = fread(buffer, 1, sizeof(buffer), file)) > 0)
    {
        char *data = buffer + numberOfBytesRead - 1;
        
        while (*data != '\n')
        {
            data--;
        }
        
        // lastNewline is our stopping point
        char *lastNewline = data;
        data = buffer;
        
        while (1)
        {
            while (data < lastNewline && *data != '\t')
            {
                data++;
            }
            
            if (data >= lastNewline)
            {
                // seek to one position after last new line character
                fseek(file, lastNewline - (buffer + numberOfBytesRead) + 1, SEEK_CUR);
                break;
            }
            
            // go past '\t'
            data++;
            
            if (listIndex >= listMaxCount)
            {
                listMaxCount *= 2;
                list = realloc(list, sizeof(float) * listMaxCount);
                if (!list)
                {
                    printf("realloc list on %s failed\n", filepath);
                    fclose(file);
                    return NULL;
                }
            }
            
            list[listIndex] = atof(data);
            listIndex++;
            
            while (data < lastNewline && *data != '\n')
            {
                data++;
            }
            // #warning, code below only works if file has two columns, versus the line above this which will work for n >= 2 columns
            // one character from float we read + one new line + one character from 1st column of next line - we can skip at least this much
            // data += 3;
        }
    }
    
    char tempPath[PATH_MAX+1];
    strncpy(tempPath, outfilepath, strlen(outfilepath)+1);
    strcat(tempPath, "_temp");
    
    cachedFile = fopen(tempPath, "w");
    if (cachedFile)
    {
        if (isByteOrderLittleEndian)
        {
            if (fwrite(list, sizeof(float), listIndex, cachedFile) < listIndex)
            {
                printf("Failed to write file %s\n", tempPath);
                return NULL;
            }
        }
        else
        {
            assert(sizeof(float) == 4);
            
            int i;
            for (i = 0; i < listIndex; i++)
            {
                float value = list[i];
                ((uint8_t *)&value)[0] = ((uint8_t *)&list[i])[3];
                ((uint8_t *)&value)[1] = ((uint8_t *)&list[i])[2];
                ((uint8_t *)&value)[2] = ((uint8_t *)&list[i])[1];
                ((uint8_t *)&value)[3] = ((uint8_t *)&list[i])[0];
                
                if (fwrite(&value, sizeof(float), 1, cachedFile) < 1)
                {
                    printf("Failed to write data to file %s\n", tempPath);
                    return NULL;
                }
            }
        }
        
        fclose(cachedFile);
        cachedFile = NULL;
        
        if (rename(tempPath, outfilepath) < 0)
        {
            printf("Failed to move %s to %s\n", tempPath, outfilepath);
            return NULL;
        }
    }
    else
    {
        printf("Writing out cached file to %s failed\n", outfilepath);
    }
    
    return Py_BuildValue("");
}

/*
static PyObject *buildList(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;
    
    freeList(self);
    
    const char *filepath = NULL;
    if (PyArg_ParseTuple(args, "s", &filepath))
    {
        FILE *file = fopen(filepath, "r");
        if (file)
        {
            unsigned int listMaxCount = 1000000000;
            
            list = malloc(sizeof(float) * listMaxCount);
            if (!list)
            {
                printf("Could not allocate enough memory.\n");
                return result;
            }
            
            // start at index 1 to match line numbers
            list[0] = 0.0;
            unsigned int listIndex = 1;
            
            char buffer[256];
            while (fgets(buffer, sizeof(buffer), file))
            {
                char *str = buffer;
                while (*str != '\t')
                {
                    str++;
                }
                str++;
                
                if (listIndex >= listMaxCount)
                {
                    listMaxCount *= 2;
                    list = realloc(list, sizeof(float) * listMaxCount);
                }
                
                list[listIndex] = atof(str);
                listIndex++;
            }
            
            fclose(file);
            // have to return something
            result = Py_BuildValue("");
        }
    }
    
    return result;
}
*/


static PyObject *floatInList(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;
    
    unsigned int index = 0;
    if (PyArg_ParseTuple(args, "I", &index))
    {
        if (cachedFile)
        {
            if (fseek(cachedFile, index * sizeof(float), SEEK_SET) < 0)
            {
                printf("Failed to seek in cached file.\n");
                return NULL;
            }
            float value = 0.0;
            if (fread(&value, sizeof(float), 1, cachedFile) < 1)
            {
                printf("Failed to read value from cached file.\n");
                return NULL;
            }
            if (!isByteOrderLittleEndian)
            {
                float value2 = value;
                ((uint8_t *)&value)[0] = ((uint8_t *)&value2)[3];
                ((uint8_t *)&value)[1] = ((uint8_t *)&value2)[2];
                ((uint8_t *)&value)[2] = ((uint8_t *)&value2)[1];
                ((uint8_t *)&value)[3] = ((uint8_t *)&value2)[0];
            }
            
            result = Py_BuildValue("f", value);
        }
        else
        {
            result = Py_BuildValue("f", list[index]);
        }
    }
    
    return result;
}

PyMethodDef gerprate_methods[] =
{
    {"buildList", (PyCFunction)buildList, METH_VARARGS, NULL},
    {"floatInList", (PyCFunction)floatInList, METH_VARARGS, NULL},
    {"freeMemory", (PyCFunction)freeMemory, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC initgerprate(void)
{
    Py_InitModule("gerprate", gerprate_methods);
    isByteOrderLittleEndian = isLittleEndian();
}
