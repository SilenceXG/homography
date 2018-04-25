import serial
port = 'COM3'
with serial.Serial('port', 9600, timeout = 0) as ser:
    #0xffffffff signifies start
    while(mode = ser.read(4) != 0xffffffff):
        pass
    # read 8 feature points
    # Mode 0: capture mode. Store Feature points as target
    # Mode 1: Navigation mode. Store Feature points as curr
    if mode == 0:
        f = open('target.txt', 'w+')
        for i in range(24):
            f.write(ser.read(4))
    elif mode == 1:
        #run exe
        pass

    
