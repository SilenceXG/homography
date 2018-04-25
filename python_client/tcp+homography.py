import numpy as np
import cv2
from matplotlib import pyplot as plt
import socket
from PIL import Image

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

cvIm = np.array(im)

# Homography Estimation
MIN_MATCH_COUNT = 10

img1 = cv2.imread('target.bmp',0)       # targetImage
img2 = np.array(im)

# Initiate SIFT detector
sift = cv2.xfeatures2d.SIFT_create()

# find the keypoints and descriptors with SIFT
kp1, des1 = sift.detectAndCompute(img1,None)
kp2, des2 = sift.detectAndCompute(img2,None)

FLANN_INDEX_KDTREE = 0
index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
search_params = dict(checks = 50)

flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1,des2,k=2)

# store all the good matches as per Lowe's ratio test.
good = []
for m,n in matches:
    if m.distance <  0.7 * n.distance:
        good.append(m)

'''
if len(good)>MIN_MATCH_COUNT:
    src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
    dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)

    M, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC,5.0)
    matchesMask = mask.ravel().tolist()

    h,w = img1.shape
    pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
    dst = cv2.perspectiveTransform(pts,M)

    img2 = cv2.polylines(img2,[np.int32(dst)],True,255,3, cv2.LINE_AA)
'''
if len(good)>MIN_MATCH_COUNT:
    src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
    dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)
    
    M, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC,5.0)
    matchesMask = mask.ravel().tolist()
    #M1, mask = cv2.findHomography(src_pts[0:(int)(len(src_pts)/2)], dst_pts[0:(int)(len(dst_pts)/2)], cv2.RANSAC,5.0)
    #M2, mask = cv2.findHomography(src_pts[(int)(len(src_pts)/2) + 1:len(src_pts)], dst_pts[(int)(len(dst_pts)/2) + 1:len(dst_pts)], cv2.LMEDS)
    M1, mask = cv2.findHomography(src_pts[1::2], dst_pts[1::2], cv2.RANSAC,5.0)
    M2, mask = cv2.findHomography(src_pts[0::2], dst_pts[0::2], cv2.RANSAC,5.0)
    h,w = img1.shape
    pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
    dst = cv2.perspectiveTransform(pts,M)

    img2 = cv2.polylines(img2,[np.int32(dst)],True,255,3, cv2.LINE_AA)
    
else:
    print ("Not enough matches are found - %d/%d" % (len(good),MIN_MATCH_COUNT))
    #translation_vector = [-1, -1, -1];
    s.close()
    matchesMask = None
 
draw_params = dict(matchColor = (0,255,0), # draw matches in green color
                   singlePointColor = None,
                   matchesMask = matchesMask, # draw only inliers
                   flags = 2)

img3 = cv2.drawMatches(img1,kp1,img2,kp2,good,None,**draw_params)

plt.imshow(img3, 'gray'),plt.show()


#homography decompostition
cameraMat = np.array([[1.0132224982983303e+03, 0.0, 320.0],[0.0, 1.0132224982983303e+03, 240],[0.0, 0.0, 1.0]])
_,Rs, Ts, Ns = cv2.decomposeHomographyMat(M1, cameraMat)
_,Rs2, Ts2, Ns2 = cv2.decomposeHomographyMat(M2, cameraMat)
_,Rs3, Ts3, Ns3 = cv2.decomposeHomographyMat(M, cameraMat)
ang1 = np.arctan(Ts[0][1]/Ts[0][0])
ang2 = np.arctan(Ts[1][1]/Ts[1][0])
ang3 = np.arctan(Ts[2][1]/Ts[2][0])
ang4 = np.arctan(Ts[3][1]/Ts[3][0])

ang12 = np.arctan(Ts2[0][1]/Ts2[0][0])
ang22 = np.arctan(Ts2[1][1]/Ts2[1][0])
ang32 = np.arctan(Ts2[2][1]/Ts2[2][0])
ang42 = np.arctan(Ts2[3][1]/Ts2[3][0])

ang13 = np.arctan(Ts3[0][1]/Ts3[0][0])
ang23 = np.arctan(Ts3[1][1]/Ts3[1][0])
ang33 = np.arctan(Ts3[2][1]/Ts3[2][0])
ang43 = np.arctan(Ts3[3][1]/Ts3[3][0])
print(ang1, ang2, ang3, ang4)
print(ang12, ang22, ang32, ang42)
print(ang13, ang23, ang33, ang43)
#TODO: Use normal vector to decide final two vectors
'''
Ts1 = Ts[0]
Ts2 = Ts[1]
ang1 = np.arctan(Ts1[1] / Ts1[0])
ang2 = np.arctan(Ts2[1] / Ts2[0])
'''
s.close()

