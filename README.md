
### rattlegram-cli

Basic C++ cli version of the Rattlegram app that uses wav files as input/output of the audio signals.

Encode usage:
```
./encode MESSAGE CALLSIGN [NOISE_SYMBOLS] [CARRIER_FREQUENCY] [RATE] [BITS] [CHANNEL] [PSK] [FILE]
```

Decode usage:
```
./decode FILE CHANNEL [PSK]
```
#### Dependencies
AudioFile https://github.com/adamstark/AudioFile
