#!/usr/bin/python

import Adafruit_DHT
import RPi.GPIO as GPIO
import ADC0832
import time

sensor = Adafruit_DHT.DHT11

# Use GPIO17 for DHT11
dht_pin = 17
LedPin = 18

Led_status = 0

def setup():
  GPIO.setmode(GPIO.BCM)       # Numbers GPIOs by physical location
  GPIO.setup(LedPin, GPIO.OUT)   # Set LedPin's mode is output

def ledOn():
  GPIO.output(LedPin, GPIO.HIGH)

def ledOff():
  GPIO.output(LedPin, GPIO.LOW)

def readTempHumidity():
  print "Reading humidate"
  humidity, temperature = Adafruit_DHT.read_retry(sensor, dht_pin)

  if humidity is not None and temperature is not None:
    print 'Temp={0:0.1f}*C  Humidity={1:0.1f}%'.format(temperature, humidity)
  else:
    print 'No reading from sensor.'

def readPhotoSens(chn=0):
  print "Reading light from channel: %d" % chn
  ADC0832.setup()
  photores = ADC0832.getResult(chn) - 80
  if photores < 0:
    photores = 0
  if photores > 100:
    photores = 100
  print 'Relative light = %d' % photores
  time.sleep(0.2)

def destroy():
  GPIO.output(LedPin, GPIO.LOW)     # led off
  GPIO.cleanup()                     # Release resource

if __name__ == '__main__':     # Program start from here
  try:
    setup()

    ledOn()
    readTempHumidity()

    readPhotoSens(0)
    ledOff()

    destroy()
  except KeyboardInterrupt:
    destroy()
