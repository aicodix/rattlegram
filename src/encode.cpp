#include "encoder.hh"
#include <cstdint>
#include <memory_resource>

static EncoderInterface *encoder;
//static DecoderInterface *decoder;
typedef DSP::Complex<float> cmplx;
int RATE = 8000;
int CHANNEL = 1;
int DEFAULT_NOISE_SYMBOLS = 6;
const char *CALL_SIGN = "ANONYMOUS";
const char *MESSAGE = "TEST";
int CARRIER_FREQUENCY = 1500;
const char *OUTPUT_FILE = "test.wav";
int OUTPUT_RATE = 8000;
int OUTPUT_BITS = 16;

int main(int argc, char **argv) {
    DSP::WriteWAV<float> output_file(OUTPUT_FILE, OUTPUT_RATE, OUTPUT_BITS, CHANNEL);
    encoder = new(std::nothrow) Encoder<8000>(&output_file);
    cmplx buffer[10*1440];
    encoder->configure(reinterpret_cast<const uint8_t*>(&MESSAGE[0]), 
        reinterpret_cast<const int8_t*>(&CALL_SIGN[0]), 
        CARRIER_FREQUENCY, 
        DEFAULT_NOISE_SYMBOLS, 
        false);
    for (int i = 0; i < 5; ++i) {
        encoder->produce_write(2);
	}
    return 0;
}