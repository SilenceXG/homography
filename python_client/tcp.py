import socket
from PIL import Image
import numpy as np
import cv2

TCP_IP = '169.254.10.123'
TCP_PORT = 36000
TOTAL_SIZE = 307200
BUFFER_SIZE = 4096
MESSAGE = "capture_frame"

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((TCP_IP, TCP_PORT))
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
#im.show()

cvIm = np.array(im)
#cv2.imshow('Camera', cvIm)

gray = np.float32(cvIm)
dst = cv2.cornerHarris(gray, 2, 3, 0.04)
dst = cv2.dilate(dst, None)

cvIm[dst > 0.0001*dst.max()] = 0xff
cv2.imshow('dst', cvIm)

if cv2.waitkey(0) & 0xff == 27:
    cv2.destroyAllWindows


'''
ret, dst = cv2.threshold(dst, 0.01*dst.max(), 255, 0)
dst = np.uint8(dst)

ret, labels, stats, centroids = cv2.connectedComponentsWithStats(dst)
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.001)
corners = cv2.cornerSubPix(gray, np.float32(centroids), (5,5), (-1, -1), criteria)
res = np.hstack(centroids, corners)
res = np.int0(res)
'''

cv2.imshow('dst', cvIm)

if cv2.waitKey(0) & 0xff == 27:
    cv2.destroyAllWindows



s.close()
