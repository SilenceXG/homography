import socket
from PIL import Image

TCP_IP = '169.254.10.123'
TCP_PORT = 36000
TOTAL_SIZE = 307200
BUFFER_SIZE = 4096
MESSAGE = "capture_frame"

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((TCP_IP,

           TCP_PORT))
byt = MESSAGE.encode()
f = open('image.bin', 'wb')
data = b''
s.send(byt)
count = 0
while True:
    part = s.recv(BUFFER_SIZE)
    data += part
    count += 1
    if not part:
        print('something')
    if len(data) >= TOTAL_SIZE:
        #print(len(part))
        break
f.write(data)
print('Write finish, count ' + str(count))
f.close()

im = Image.frombytes('L', (640,480), data)
im.show()
im.save('target.bmp')
s.close()
