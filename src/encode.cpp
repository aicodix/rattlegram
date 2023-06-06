#include "encoder.hh"
#include <cstdint>

static EncoderInterface *encoder;
const int RATE = 8000;
const int CHANNEL = 3;
const int DEFAULT_NOISE_SYMBOLS = 1;
const int CARRIER_FREQUENCY = 1500;
const int OUTPUT_BITS = 16;

int main(int argc, char **argv) {
    char message[192] = {0};
    std::string mesg = "TEST";
    for (int i = 0; i < mesg.length(); i++) {
        message[i] = mesg[i];
    }
    const char call_sign[] = "OWO";
    const char out_file[] = "out.wav";
    int channel_count = 1;
    if (CHANNEL != 0) {
        channel_count = 2;
    }
    DSP::WriteWAV<float> output_file(out_file, RATE, OUTPUT_BITS, channel_count);
    encoder = new(std::nothrow) Encoder<RATE>(&output_file);
    encoder->configure(reinterpret_cast<const uint8_t*>(&message[0]), 
        reinterpret_cast<const int8_t*>(&call_sign[0]), 
        CARRIER_FREQUENCY, 
        DEFAULT_NOISE_SYMBOLS, 
        false);
    for (int i = 0; i < 9 + DEFAULT_NOISE_SYMBOLS; ++i) {
        encoder->produce_write(CHANNEL);
	}

    return 0;
}