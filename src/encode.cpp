#include "encoder.hh"
#include <cstdint>

static EncoderInterface *encoder;
//static DecoderInterface *decoder;
typedef DSP::Complex<float> cmplx;
const int RATE = 8000;
const int SYMBOL_LENGHT = (1280 * RATE) / 8000;
const int GUARD_LENGTH = SYMBOL_LENGHT/8;
const int EXTENDED_LENGTH = SYMBOL_LENGHT + GUARD_LENGTH;
const int CHANNEL = 0;
const int NUM_CHANNELS = 1;
const int DEFAULT_NOISE_SYMBOLS = 1;
const int CARRIER_FREQUENCY = 1500;
const int OUTPUT_RATE = 8000;
const int OUTPUT_BITS = 16;

int main(int argc, char **argv) {
    char MESSAGE[192] = {0};
    std::string mesg = "TEST";
    for (int i = 0; i < mesg.length(); i++) {
        MESSAGE[i] = mesg[i];
    }
    const char CALL_SIGN[] = "OWO";
    const char OUTPUT_FILE[] = "test.wav";
    DSP::WriteWAV<float> output_file(OUTPUT_FILE, OUTPUT_RATE, OUTPUT_BITS, NUM_CHANNELS);
    encoder = new(std::nothrow) Encoder<8000>(&output_file);
    int16_t buffer[10*EXTENDED_LENGTH] = {0};
    encoder->configure(reinterpret_cast<const uint8_t*>(&MESSAGE[0]), 
        reinterpret_cast<const int8_t*>(&CALL_SIGN[0]), 
        CARRIER_FREQUENCY, 
        DEFAULT_NOISE_SYMBOLS, 
        false);
    for (int i = 0; i < 10; ++i) {
        //encoder->produce(&buffer[i*NUM_CHANNELS*EXTENDED_LENGTH], CHANNEL);
        encoder->produce_write(1);
	}

    return 0;
}