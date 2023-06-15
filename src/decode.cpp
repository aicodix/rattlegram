#include "decoder.hh"
#include "wav.hh"
#include <cstdint>
#include "AudioFile.h"

static DecoderInterface *decoder;

int main(int argc, char **argv) {
	char input_name[] = "out.wav";
	AudioFile<int16_t> audioFile;
	audioFile.load(input_name);
	int rate = audioFile.getSampleRate();
	int channel_count = audioFile.getNumChannels();
	int file_length = audioFile.getNumSamplesPerChannel() * channel_count;
	int symbol_length = (1280 * rate) / 8000;
	int guard_length = symbol_length / 8;
	int extended_length = symbol_length + guard_length;
	int channel = 0;
	int record_count = rate/50;

	int16_t file[file_length + record_count];
	for (int i = 0; i < file_length + record_count; i++) {
		file[i] = 0;
	}
	for (int i = 0; i < file_length/channel_count; i++) {
		for (int c = 0; c < channel_count; c++) {
			file[i * channel_count + c] = audioFile.samples[c][i];
		}
	}

    switch (rate) {
		case 8000:
			decoder = new(std::nothrow) Decoder<8000>();
			break;
		case 16000:
			decoder = new(std::nothrow) Decoder<16000>();
			break;
		case 32000:
			decoder = new(std::nothrow) Decoder<32000>();
			break;
		case 44100:
			decoder = new(std::nothrow) Decoder<44100>();
			break;
		case 48000:
			decoder = new(std::nothrow) Decoder<48000>();
			break;
		default:
            std::cerr << "Unsupported sample rate." ;
            std::cerr << "Supported rates: 8000/16000/32000/44100/48000."<< std::endl;
            return 1;
	}

	for (int i = 0; i * record_count < file_length; i++) {
		if (decoder->feed(&file[i*record_count*channel_count], record_count*channel_count, channel)) {
			std::cout << i << std::endl;
			int status = decoder->process();
			std::cout << "Status: " << status << std::endl;	
			float cfo = -1.0;
			int32_t mode = -1;
			uint8_t call_sign[192] = {0}; 
			uint8_t payload[192] = {0};
			std::string call;
			std::string pay;

			switch (status) {
				case 0:
					std::cout << "OKAY" << std::endl;
					break;
				case 1:
					std::cout << "PREAMBLE FAIL" << std::endl;
					break;
				case 2:
					decoder->staged(&cfo, &mode, call_sign);
					call = (reinterpret_cast<char*>(call_sign));
					std::cout << "SYNC:"
							<< "  CFO "
							<< cfo
							<< "  mode: "
							<< mode
							<< "  call sign: "
							<< call
							<< std::endl;
					break;
				case 3:
					decoder->fetch(payload);
					pay = (reinterpret_cast<char*>(payload));
					std::cout << "DONE:  payload: " 
							<< pay
							<< std::endl;
					break;
				case 4:
					std::cout << "HEAP ERROR" << std::endl;
					break;
				case 5:
					decoder->staged(&cfo, &mode, call_sign);
					call = (reinterpret_cast<char*>(call_sign));
					std::cout << "NOPE:"
							<< "  CFO "
							<< cfo
							<< "  mode: "
							<< mode
							<< "  call sign: "
							<< call
							<< std::endl;
					break;
				case 6: 
					decoder->staged(&cfo, &mode, call_sign);
					call = (reinterpret_cast<char*>(call_sign));
					std::cout << "PING:"
							<< "  CFO "
							<< cfo
							<< "  mode: "
							<< mode
							<< "  call sign: "
							<< call
							<< std::endl;
					break;
				default:
					std::cout << "UNDEFINED BEHAVIOUR" << std::endl;
					return 1;

			}
		}
	}
	uint8_t payload[192] = {0};
	std::string payload_str;
	decoder->fetch(payload);
	payload_str = (reinterpret_cast<char*>(payload));
	std::cout << "DONE:  payload: " 
			<< payload_str
			<< std::endl;
    return 0;
}