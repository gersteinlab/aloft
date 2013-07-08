#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

// Returns 1 if the machine is little endian, otherwise 0 if big endian
static char isLittleEndian(void)
{
    int testInteger = 1;
    return (*((unsigned char *)(&testInteger)) == 1);
}

static void writeFloat(FILE *file, char isLittleEndianByteOrder, float data)
{
    if (isLittleEndianByteOrder)
    {
        if (fwrite(&data, sizeof(data), 1, file) < 1)
        {
            printf("Failed to write float to file..\n");
            exit(1);
        }
    }
    else
    {
        float value = data;
        ((uint8_t *)&value)[0] = ((uint8_t *)&data)[3];
        ((uint8_t *)&value)[1] = ((uint8_t *)&data)[2];
        ((uint8_t *)&value)[2] = ((uint8_t *)&data)[1];
        ((uint8_t *)&value)[3] = ((uint8_t *)&data)[0];
        
        if (fwrite(&value, sizeof(float), 1, file) < 1)
        {
            printf("Failed to write little endian data to file..\n");
            exit(1);
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        printf("Usage: %s <gerp_input_file> <cache_output_file>\n", argv[0]);
        exit(1);
    }
    
    const char *filepath = argv[1];
    const char *outfilepath = argv[2];

    char isLittleEndianByteOrder = isLittleEndian();
    
    FILE *file = fopen(filepath, "r");
    if (!file)
    {
        printf("Failed to open file %s\n", filepath);
        exit(1);
    }

    char tempPath[PATH_MAX+1];
    strncpy(tempPath, outfilepath, strlen(outfilepath)+1);
    strcat(tempPath, "_temp");

    FILE *cachedFile = fopen(tempPath, "w");
    if (!cachedFile)
    {
        printf("Failed to open file for writing: %s\n", tempPath);
        exit(1);
    }

    // start at index 1 to match line numbers
    writeFloat(cachedFile, isLittleEndianByteOrder, 0.0);
    
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
            
            writeFloat(cachedFile, isLittleEndianByteOrder, atof(data));
            
            while (data < lastNewline && *data != '\n')
            {
                data++;
            }
            // #warning, code below only works if file has two columns, versus the line above this which will work for n >= 2 columns
            // one character from float we read + one new line + one character from 1st column of next line - we can skip at least this much
            // data += 3;
        }
    }
    
    if (rename(tempPath, outfilepath) < 0)
    {
        printf("Failed to move %s to %s\n", tempPath, outfilepath);
        exit(1);
    }

    fclose(file);
    fclose(cachedFile);
    
    return 0;
}
