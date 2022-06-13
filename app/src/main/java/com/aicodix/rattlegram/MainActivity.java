/*
Rattlegram

Copyright 2022 Ahmet Inan <inan@aicodix.de>
*/

package com.aicodix.rattlegram;

import androidx.annotation.NonNull;
import androidx.appcompat.app.AlertDialog;
import androidx.appcompat.app.AppCompatActivity;
import androidx.appcompat.app.AppCompatDelegate;
import androidx.core.app.ActivityCompat;
import androidx.core.content.ContextCompat;

import android.Manifest;
import android.content.Context;
import android.content.SharedPreferences;
import android.content.pm.PackageManager;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioRecord;
import android.media.AudioTrack;
import android.media.MediaRecorder;
import android.os.Build;
import android.os.Bundle;
import android.text.InputType;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ArrayAdapter;
import android.widget.EditText;
import android.widget.NumberPicker;

import com.aicodix.rattlegram.databinding.ActivityMainBinding;

import java.nio.charset.StandardCharsets;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Locale;

public class MainActivity extends AppCompatActivity {

	// Used to load the 'rattlegram' library on application startup.
	static {
		System.loadLibrary("rattlegram");
	}

	private final int permissionID = 1;
	private final int operationMode = 14;
	private final int audioFormat = AudioFormat.ENCODING_PCM_16BIT;
	private final int sampleSize = 2;
	private AudioRecord audioRecord;
	private AudioTrack audioTrack;
	private boolean fancyHeader;
	private boolean repeaterMode;
	private int noiseSymbols;
	private int recordRate;
	private int outputRate;
	private int recordChannel;
	private int outputChannel;
	private int audioSource;
	private int carrierFrequency;
	private int bufferLength;
	private short[] recordBuffer;
	private short[] outputBuffer;
	private Menu menu;
	private byte[] payload;
	private ArrayAdapter<String> messages;
	private float[] cachedCFO;
	private int[] cachedMode;
	private byte[] cachedCall;
	private String callSign;
	private String draft;

	private native boolean createEncoder(int sampleRate);

	private native void configureEncoder(byte[] payload, byte[] callSign, int carrierFrequency, int noiseSymbols, boolean fancyHeader);

	private native boolean produceEncoder(short[] audioBuffer, int channelSelect);

	private native void destroyEncoder();

	private final AudioTrack.OnPlaybackPositionUpdateListener outputListener = new AudioTrack.OnPlaybackPositionUpdateListener() {
		@Override
		public void onMarkerReached(AudioTrack ignore) {

		}

		@Override
		public void onPeriodicNotification(AudioTrack audioTrack) {
			if (produceEncoder(outputBuffer, outputChannel)) {
				audioTrack.write(outputBuffer, 0, outputBuffer.length);
			} else {
				audioTrack.stop();
				startListening();
			}
		}
	};

	private void initAudioTrack() {
		if (audioTrack != null) {
			boolean rateChanged = audioTrack.getSampleRate() != outputRate;
			boolean channelChanged = audioTrack.getChannelCount() != (outputChannel == 0 ? 1 : 2);
			if (!rateChanged && !channelChanged)
				return;
			audioTrack.stop();
			audioTrack.release();
		}
		int channelConfig = AudioFormat.CHANNEL_OUT_MONO;
		int channelCount = 1;
		if (outputChannel != 0) {
			channelCount = 2;
			channelConfig = AudioFormat.CHANNEL_OUT_STEREO;
		}
		bufferLength = 2 * Integer.highestOneBit(outputRate) * channelCount;
		int bufferSize = sampleSize * bufferLength;
		int symbolLength = (1280 * outputRate) / 8000;
		int guardLength = symbolLength / 8;
		int extendedLength = symbolLength + guardLength;
		audioTrack = new AudioTrack(AudioManager.STREAM_MUSIC, outputRate, channelConfig, audioFormat, bufferSize, AudioTrack.MODE_STREAM);
		outputBuffer = new short[extendedLength * channelCount];
		audioTrack.setPlaybackPositionUpdateListener(outputListener);
		audioTrack.setPositionNotificationPeriod(extendedLength);
		if (!createEncoder(outputRate))
			addMessage(getString(R.string.encoder_status), getString(R.string.heap_error));
	}

	private native int processDecoder(short[] audioBuffer, int channelSelect);

	private native void cachedDecoder(float[] carrierFrequencyOffset, int[] operationMode, byte[] callSign);

	private native boolean fetchDecoder(byte[] payload);

	private native boolean createDecoder(int sampleRate);

	private native void destroyDecoder();

	private final AudioRecord.OnRecordPositionUpdateListener recordListener = new AudioRecord.OnRecordPositionUpdateListener() {
		@Override
		public void onMarkerReached(AudioRecord ignore) {

		}

		@Override
		public void onPeriodicNotification(AudioRecord audioRecord) {
			audioRecord.read(recordBuffer, 0, recordBuffer.length);
			int status = processDecoder(recordBuffer, recordChannel);
			final int STATUS_OKAY = 0;
			final int STATUS_FAIL = 1;
			final int STATUS_SYNC = 2;
			final int STATUS_DONE = 3;
			final int STATUS_HEAP = 4;
			final int STATUS_NOPE = 5;
			switch (status) {
				case STATUS_OKAY:
					break;
				case STATUS_FAIL:
					addMessage(getString(R.string.decoder_status), getString(R.string.preamble_fail));
					break;
				case STATUS_NOPE:
					cachedDecoder(cachedCFO, cachedMode, cachedCall);
					addMessage(fromTitle(), getString(R.string.preamble_nope));
					break;
				case STATUS_HEAP:
					addMessage(getString(R.string.decoder_status), getString(R.string.heap_error));
					break;
				case STATUS_SYNC:
					cachedDecoder(cachedCFO, cachedMode, cachedCall);
					break;
				case STATUS_DONE:
					if (fetchDecoder(payload)) {
						addMessage(fromTitle(), new String(payload).trim());
						if (repeaterMode)
							repeatMessage();
					} else {
						addMessage(fromTitle(), getString(R.string.decoding_failed));
					}
					break;
			}
		}
	};

	private byte[] callTerm() {
		return Arrays.copyOf(callSign.getBytes(StandardCharsets.US_ASCII), callSign.length() + 1);
	}

	private String currentTime() {
		return new SimpleDateFormat("HH:mm:ss", Locale.US).format(new Date());
	}

	private String fromTitle() {
		return getString(R.string.from_title, new String(cachedCall).trim(), cachedMode[0], cachedCFO[0]);
	}

	private String sentTitle() {
		return getString(R.string.sent_title, callSign.trim(), operationMode, carrierFrequency);
	}

	private String repeatTitle() {
		return getString(R.string.repeat_title, new String(cachedCall).trim(), operationMode, carrierFrequency);
	}

	private void addMessage(String title, String mesg) {
		addString(getString(R.string.title_message, currentTime(), title, mesg));
	}

	private void addString(String str) {
		int count = 100;
		if (messages.getCount() >= count)
			messages.remove(messages.getItem(count - 1));
		messages.insert(str, 0);
	}

	private void startListening() {
		if (audioRecord != null) {
			audioRecord.startRecording();
			if (audioRecord.getRecordingState() == AudioRecord.RECORDSTATE_RECORDING) {
				audioRecord.read(recordBuffer, 0, recordBuffer.length);
			} else {
				addMessage(getString(R.string.recorder_status), getString(R.string.audio_recording_error));
			}
		}
	}

	private void stopListening() {
		if (audioRecord != null)
			audioRecord.stop();
	}

	private void initAudioRecord(boolean restart) {
		if (audioRecord != null) {
			boolean rateChanged = audioRecord.getSampleRate() != recordRate;
			boolean channelChanged = audioRecord.getChannelCount() != (recordChannel == 0 ? 1 : 2);
			boolean sourceChanged = audioRecord.getAudioSource() != audioSource;
			if (!rateChanged && !channelChanged && !sourceChanged)
				return;
			stopListening();
			audioRecord.release();
			audioRecord = null;
		}
		int channelConfig = AudioFormat.CHANNEL_IN_MONO;
		int channelCount = 1;
		if (recordChannel != 0) {
			channelCount = 2;
			channelConfig = AudioFormat.CHANNEL_IN_STEREO;
		}
		int frameSize = sampleSize * channelCount;
		int bufferSize = 2 * Integer.highestOneBit(3 * recordRate) * frameSize;
		int symbolLength = (1280 * recordRate) / 8000;
		int guardLength = symbolLength / 8;
		int extendedLength = symbolLength + guardLength;
		try {
			AudioRecord testAudioRecord = new AudioRecord(audioSource, recordRate, channelConfig, audioFormat, bufferSize);
			if (testAudioRecord.getState() == AudioRecord.STATE_INITIALIZED) {
				if (createDecoder(recordRate)) {
					audioRecord = testAudioRecord;
					recordBuffer = new short[extendedLength * channelCount];
					audioRecord.setRecordPositionUpdateListener(recordListener);
					audioRecord.setPositionNotificationPeriod(extendedLength);
					if (restart)
						startListening();
				} else {
					testAudioRecord.release();
					addMessage(getString(R.string.encoder_status), getString(R.string.heap_error));
				}
			} else {
				testAudioRecord.release();
				addMessage(getString(R.string.recorder_status), getString(R.string.audio_init_failed));
			}
		} catch (IllegalArgumentException e) {
			addMessage(getString(R.string.recorder_status), getString(R.string.audio_setup_failed));
		} catch (SecurityException e) {
			addMessage(getString(R.string.recorder_status), getString(R.string.audio_permission_denied));
		}
	}

	private void setRecordRate(int newSampleRate) {
		if (recordRate == newSampleRate)
			return;
		recordRate = newSampleRate;
		updateRecordRateMenu();
		initAudioRecord(true);
	}

	private void setRecordChannel(int newChannelSelect) {
		if (recordChannel == newChannelSelect)
			return;
		recordChannel = newChannelSelect;
		updateRecordChannelMenu();
		initAudioRecord(true);
	}

	private void setAudioSource(int newAudioSource) {
		if (audioSource == newAudioSource)
			return;
		audioSource = newAudioSource;
		updateAudioSourceMenu();
		initAudioRecord(true);
	}

	@Override
	public void onRequestPermissionsResult(int requestCode, @NonNull String[] permissions, @NonNull int[] grantResults) {
		super.onRequestPermissionsResult(requestCode, permissions, grantResults);
		if (requestCode != permissionID)
			return;
		for (int i = 0; i < permissions.length; ++i)
			if (permissions[i].equals(Manifest.permission.RECORD_AUDIO) && grantResults[i] == PackageManager.PERMISSION_GRANTED)
				initAudioRecord(false);
	}

	@Override
	protected void onSaveInstanceState(@NonNull Bundle state) {
		state.putInt("nightMode", AppCompatDelegate.getDefaultNightMode());
		state.putInt("outputRate", outputRate);
		state.putInt("outputChannel", outputChannel);
		state.putInt("recordRate", recordRate);
		state.putInt("recordChannel", recordChannel);
		state.putInt("audioSource", audioSource);
		state.putInt("carrierFrequency", carrierFrequency);
		state.putInt("noiseSymbols", noiseSymbols);
		state.putString("callSign", callSign);
		state.putBoolean("fancyHeader", fancyHeader);
		state.putBoolean("repeaterMode", repeaterMode);
		super.onSaveInstanceState(state);
	}

	private void storeSettings() {
		SharedPreferences pref = getPreferences(Context.MODE_PRIVATE);
		SharedPreferences.Editor edit = pref.edit();
		edit.putInt("nightMode", AppCompatDelegate.getDefaultNightMode());
		edit.putInt("outputRate", outputRate);
		edit.putInt("outputChannel", outputChannel);
		edit.putInt("recordRate", recordRate);
		edit.putInt("recordChannel", recordChannel);
		edit.putInt("audioSource", audioSource);
		edit.putInt("carrierFrequency", carrierFrequency);
		edit.putInt("noiseSymbols", noiseSymbols);
		edit.putString("callSign", callSign);
		edit.putBoolean("fancyHeader", fancyHeader);
		edit.putBoolean("repeaterMode", repeaterMode);
		edit.apply();
	}

	@Override
	protected void onCreate(Bundle state) {
		final int defaultSampleRate = 8000;
		final int defaultChannelSelect = 0;
		final int defaultAudioSource = MediaRecorder.AudioSource.DEFAULT;
		final int defaultCarrierFrequency = 1450;
		final int defaultNoiseSymbols = 6;
		final String defaultCallSign = "ANONYMOUS";
		final boolean defaultFancyHeader = false;
		final boolean defaultRepeaterMode = false;
		if (state == null) {
			SharedPreferences pref = getPreferences(Context.MODE_PRIVATE);
			AppCompatDelegate.setDefaultNightMode(pref.getInt("nightMode", AppCompatDelegate.getDefaultNightMode()));
			outputRate = pref.getInt("outputRate", defaultSampleRate);
			outputChannel = pref.getInt("outputChannel", defaultChannelSelect);
			recordRate = pref.getInt("recordRate", defaultSampleRate);
			recordChannel = pref.getInt("recordChannel", defaultChannelSelect);
			audioSource = pref.getInt("audioSource", defaultAudioSource);
			carrierFrequency = pref.getInt("carrierFrequency", defaultCarrierFrequency);
			noiseSymbols = pref.getInt("noiseSymbols", defaultNoiseSymbols);
			callSign = pref.getString("callSign", defaultCallSign);
			fancyHeader = pref.getBoolean("fancyHeader", defaultFancyHeader);
			repeaterMode = pref.getBoolean("repeaterMode", defaultRepeaterMode);
		} else {
			AppCompatDelegate.setDefaultNightMode(state.getInt("nightMode", AppCompatDelegate.getDefaultNightMode()));
			outputRate = state.getInt("outputRate", defaultSampleRate);
			outputChannel = state.getInt("outputChannel", defaultChannelSelect);
			recordRate = state.getInt("recordRate", defaultSampleRate);
			recordChannel = state.getInt("recordChannel", defaultChannelSelect);
			audioSource = state.getInt("audioSource", defaultAudioSource);
			carrierFrequency = state.getInt("carrierFrequency", defaultCarrierFrequency);
			noiseSymbols = state.getInt("noiseSymbols", defaultNoiseSymbols);
			callSign = state.getString("callSign", defaultCallSign);
			fancyHeader = state.getBoolean("fancyHeader", defaultFancyHeader);
			repeaterMode = state.getBoolean("repeaterMode", defaultRepeaterMode);
		}
		super.onCreate(state);
		ActivityMainBinding binding = ActivityMainBinding.inflate(getLayoutInflater());
		setContentView(binding.getRoot());
		cachedCFO = new float[1];
		cachedMode = new int[1];
		cachedCall = new byte[10];
		payload = new byte[170];
		messages = new ArrayAdapter<>(this, android.R.layout.simple_list_item_1);
		binding.messages.setAdapter(messages);
		initAudioTrack();

		List<String> permissions = new ArrayList<>();
		if (ContextCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO) != PackageManager.PERMISSION_GRANTED) {
			permissions.add(Manifest.permission.RECORD_AUDIO);
			addMessage(getString(R.string.recorder_status), getString(R.string.audio_permission_denied));
		} else {
			initAudioRecord(false);
		}
		if (!permissions.isEmpty())
			ActivityCompat.requestPermissions(this, permissions.toArray(new String[0]), permissionID);
	}

	private void setNoiseSymbols(int newNoiseSymbols) {
		if (noiseSymbols == newNoiseSymbols)
			return;
		noiseSymbols = newNoiseSymbols;
		updateNoiseSymbolsMenu();
	}

	private void updateNoiseSymbolsMenu() {
		switch (noiseSymbols) {
			case 0:
				menu.findItem(R.id.action_disable_noise).setChecked(true);
				break;
			case 1:
				menu.findItem(R.id.action_set_noise_quarter_second).setChecked(true);
				break;
			case 3:
				menu.findItem(R.id.action_set_noise_half_second).setChecked(true);
				break;
			case 6:
				menu.findItem(R.id.action_set_noise_one_second).setChecked(true);
				break;
			case 11:
				menu.findItem(R.id.action_set_noise_two_seconds).setChecked(true);
				break;
			case 22:
				menu.findItem(R.id.action_set_noise_four_seconds).setChecked(true);
				break;
		}
	}

	private void setFancyHeader(boolean newFancyHeader) {
		if (fancyHeader == newFancyHeader)
			return;
		fancyHeader = newFancyHeader;
		updateFancyHeaderMenu();
	}

	private void updateFancyHeaderMenu() {
		if (fancyHeader)
			menu.findItem(R.id.action_enable_fancy_header).setChecked(true);
		else
			menu.findItem(R.id.action_disable_fancy_header).setChecked(true);
	}

	private void setRepeaterMode(boolean newRepeaterMode) {
		if (repeaterMode == newRepeaterMode)
			return;
		repeaterMode = newRepeaterMode;
		updateRepeaterModeMenu();
	}

	private void updateRepeaterModeMenu() {
		if (repeaterMode)
			menu.findItem(R.id.action_enable_repeater_mode).setChecked(true);
		else
			menu.findItem(R.id.action_disable_repeater_mode).setChecked(true);
	}

	private void setOutputRate(int newSampleRate) {
		if (audioTrack.getPlayState() == AudioTrack.PLAYSTATE_PLAYING)
			return;
		if (outputRate == newSampleRate)
			return;
		outputRate = newSampleRate;
		updateOutputRateMenu();
		initAudioTrack();
	}

	private void updateOutputRateMenu() {
		switch (outputRate) {
			case 8000:
				menu.findItem(R.id.action_set_output_rate_8000).setChecked(true);
				break;
			case 16000:
				menu.findItem(R.id.action_set_output_rate_16000).setChecked(true);
				break;
			case 32000:
				menu.findItem(R.id.action_set_output_rate_32000).setChecked(true);
				break;
			case 44100:
				menu.findItem(R.id.action_set_output_rate_44100).setChecked(true);
				break;
			case 48000:
				menu.findItem(R.id.action_set_output_rate_48000).setChecked(true);
				break;
		}
	}

	private void updateRecordRateMenu() {
		switch (recordRate) {
			case 8000:
				menu.findItem(R.id.action_set_record_rate_8000).setChecked(true);
				break;
			case 16000:
				menu.findItem(R.id.action_set_record_rate_16000).setChecked(true);
				break;
			case 32000:
				menu.findItem(R.id.action_set_record_rate_32000).setChecked(true);
				break;
			case 44100:
				menu.findItem(R.id.action_set_record_rate_44100).setChecked(true);
				break;
			case 48000:
				menu.findItem(R.id.action_set_record_rate_48000).setChecked(true);
				break;
		}
	}

	private void updateRecordChannelMenu() {
		switch (recordChannel) {
			case 0:
				menu.findItem(R.id.action_set_record_channel_default).setChecked(true);
				break;
			case 1:
				menu.findItem(R.id.action_set_record_channel_first).setChecked(true);
				break;
			case 2:
				menu.findItem(R.id.action_set_record_channel_second).setChecked(true);
				break;
			case 3:
				menu.findItem(R.id.action_set_record_channel_summation).setChecked(true);
				break;
			case 4:
				menu.findItem(R.id.action_set_record_channel_analytic).setChecked(true);
				break;
		}
	}

	private void updateAudioSourceMenu() {
		switch (audioSource) {
			case MediaRecorder.AudioSource.DEFAULT:
				menu.findItem(R.id.action_set_source_default).setChecked(true);
				break;
			case MediaRecorder.AudioSource.MIC:
				menu.findItem(R.id.action_set_source_microphone).setChecked(true);
				break;
			case MediaRecorder.AudioSource.CAMCORDER:
				menu.findItem(R.id.action_set_source_camcorder).setChecked(true);
				break;
			case MediaRecorder.AudioSource.VOICE_RECOGNITION:
				menu.findItem(R.id.action_set_source_voice_recognition).setChecked(true);
				break;
			case MediaRecorder.AudioSource.UNPROCESSED:
				menu.findItem(R.id.action_set_source_unprocessed).setChecked(true);
				break;
		}
	}

	private void setOutputChannel(int newChannelSelect) {
		if (audioTrack.getPlayState() == AudioTrack.PLAYSTATE_PLAYING)
			return;
		if (outputChannel == newChannelSelect)
			return;
		outputChannel = newChannelSelect;
		updateOutputChannelMenu();
		initAudioTrack();
	}

	private void updateOutputChannelMenu() {
		switch (outputChannel) {
			case 0:
				menu.findItem(R.id.action_set_output_channel_default).setChecked(true);
				break;
			case 1:
				menu.findItem(R.id.action_set_output_channel_first).setChecked(true);
				break;
			case 2:
				menu.findItem(R.id.action_set_output_channel_second).setChecked(true);
				break;
			case 4:
				menu.findItem(R.id.action_set_output_channel_analytic).setChecked(true);
				break;
		}
	}

	@Override
	public boolean onCreateOptionsMenu(Menu menu) {
		getMenuInflater().inflate(R.menu.menu_main, menu);
		menu.findItem(R.id.action_set_source_unprocessed).setEnabled(Build.VERSION.SDK_INT >= Build.VERSION_CODES.N);
		this.menu = menu;
		updateOutputRateMenu();
		updateOutputChannelMenu();
		updateRecordRateMenu();
		updateRecordChannelMenu();
		updateAudioSourceMenu();
		updateNoiseSymbolsMenu();
		updateFancyHeaderMenu();
		updateRepeaterModeMenu();
		return true;
	}

	@Override
	public boolean onOptionsItemSelected(MenuItem item) {
		int id = item.getItemId();
		if (id == R.id.action_compose) {
			composeMessage();
			return true;
		}
		if (id == R.id.action_set_output_rate_8000) {
			setOutputRate(8000);
			return true;
		}
		if (id == R.id.action_set_output_rate_16000) {
			setOutputRate(16000);
			return true;
		}
		if (id == R.id.action_set_output_rate_32000) {
			setOutputRate(32000);
			return true;
		}
		if (id == R.id.action_set_output_rate_44100) {
			setOutputRate(44100);
			return true;
		}
		if (id == R.id.action_set_output_rate_48000) {
			setOutputRate(48000);
			return true;
		}
		if (id == R.id.action_set_output_channel_default) {
			setOutputChannel(0);
			return true;
		}
		if (id == R.id.action_set_output_channel_first) {
			setOutputChannel(1);
			return true;
		}
		if (id == R.id.action_set_output_channel_second) {
			setOutputChannel(2);
			return true;
		}
		if (id == R.id.action_set_output_channel_analytic) {
			setOutputChannel(4);
			return true;
		}
		if (id == R.id.action_set_record_rate_8000) {
			setRecordRate(8000);
			return true;
		}
		if (id == R.id.action_set_record_rate_16000) {
			setRecordRate(16000);
			return true;
		}
		if (id == R.id.action_set_record_rate_32000) {
			setRecordRate(32000);
			return true;
		}
		if (id == R.id.action_set_record_rate_44100) {
			setRecordRate(44100);
			return true;
		}
		if (id == R.id.action_set_record_rate_48000) {
			setRecordRate(48000);
			return true;
		}
		if (id == R.id.action_set_record_channel_default) {
			setRecordChannel(0);
			return true;
		}
		if (id == R.id.action_set_record_channel_first) {
			setRecordChannel(1);
			return true;
		}
		if (id == R.id.action_set_record_channel_second) {
			setRecordChannel(2);
			return true;
		}
		if (id == R.id.action_set_record_channel_summation) {
			setRecordChannel(3);
			return true;
		}
		if (id == R.id.action_set_record_channel_analytic) {
			setRecordChannel(4);
			return true;
		}
		if (id == R.id.action_set_source_default) {
			setAudioSource(MediaRecorder.AudioSource.DEFAULT);
			return true;
		}
		if (id == R.id.action_set_source_microphone) {
			setAudioSource(MediaRecorder.AudioSource.MIC);
			return true;
		}
		if (id == R.id.action_set_source_camcorder) {
			setAudioSource(MediaRecorder.AudioSource.CAMCORDER);
			return true;
		}
		if (id == R.id.action_set_source_voice_recognition) {
			setAudioSource(MediaRecorder.AudioSource.VOICE_RECOGNITION);
			return true;
		}
		if (id == R.id.action_set_source_unprocessed) {
			if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.N) {
				setAudioSource(MediaRecorder.AudioSource.UNPROCESSED);
				return true;
			}
			return false;
		}
		if (id == R.id.action_edit_call_sign) {
			editCallSign();
			return true;
		}
		if (id == R.id.action_set_carrier_frequency) {
			setCarrierFrequency();
			return true;
		}
		if (id == R.id.action_disable_noise) {
			setNoiseSymbols(0);
			return true;
		}
		if (id == R.id.action_set_noise_quarter_second) {
			setNoiseSymbols(1);
			return true;
		}
		if (id == R.id.action_set_noise_half_second) {
			setNoiseSymbols(3);
			return true;
		}
		if (id == R.id.action_set_noise_one_second) {
			setNoiseSymbols(6);
			return true;
		}
		if (id == R.id.action_set_noise_two_seconds) {
			setNoiseSymbols(11);
			return true;
		}
		if (id == R.id.action_set_noise_four_seconds) {
			setNoiseSymbols(22);
			return true;
		}
		if (id == R.id.action_enable_fancy_header) {
			setFancyHeader(true);
			return true;
		}
		if (id == R.id.action_disable_fancy_header) {
			setFancyHeader(false);
			return true;
		}
		if (id == R.id.action_enable_repeater_mode) {
			setRepeaterMode(true);
			return true;
		}
		if (id == R.id.action_disable_repeater_mode) {
			setRepeaterMode(false);
			return true;
		}
		if (id == R.id.action_enable_night_mode) {
			AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_YES);
			return true;
		}
		if (id == R.id.action_disable_night_mode) {
			AppCompatDelegate.setDefaultNightMode(AppCompatDelegate.MODE_NIGHT_NO);
			return true;
		}
		if (id == R.id.action_force_quit) {
			storeSettings();
			System.exit(0);
			return true;
		}
		if (id == R.id.action_privacy_policy) {
			showTextPage(getString(R.string.privacy_policy), getString(R.string.privacy_policy_text));
			return true;
		}
		if (id == R.id.action_about) {
			showTextPage(getString(R.string.about), getString(R.string.about_text, BuildConfig.VERSION_NAME));
			return true;
		}
		return super.onOptionsItemSelected(item);
	}

	private void composeMessage() {
		View view = getLayoutInflater().inflate(R.layout.compose_message, null);
		EditText edit = view.findViewById(R.id.message);
		edit.setText(draft);
		draft = "";
		AlertDialog.Builder builder = new AlertDialog.Builder(this, R.style.Theme_AlertDialog);
		builder.setTitle(R.string.compose_message);
		builder.setView(view);
		builder.setNeutralButton(R.string.draft, (dialogInterface, i) -> draft = edit.getText().toString());
		builder.setNegativeButton(R.string.discard, null);
		builder.setPositiveButton(R.string.transmit, (dialogInterface, i) -> transmitMessage(edit.getText().toString()));
		builder.setOnCancelListener(dialogInterface -> draft = edit.getText().toString());
		builder.show();
	}

	private void transmitMessage(String message) {
		stopListening();
		byte[] mesg = Arrays.copyOf(message.getBytes(StandardCharsets.UTF_8), payload.length);
		addMessage(sentTitle(), new String(mesg).trim());
		configureEncoder(mesg, callTerm(), carrierFrequency, noiseSymbols, fancyHeader);
		audioTrack.write(new short[bufferLength], 0, bufferLength);
		audioTrack.play();
	}

	private void repeatMessage() {
		stopListening();
		addMessage(repeatTitle(), new String(payload).trim());
		configureEncoder(payload, cachedCall, carrierFrequency, noiseSymbols, fancyHeader);
		audioTrack.write(new short[bufferLength], 0, bufferLength);
		audioTrack.play();
	}

	private void setInputType(ViewGroup np, int it) {
		int count = np.getChildCount();
		for (int i = 0; i < count; i++) {
			final View child = np.getChildAt(i);
			if (child instanceof ViewGroup) {
				setInputType((ViewGroup) child, it);
			} else if (child instanceof EditText) {
				EditText et = (EditText) child;
				et.setInputType(it);
				break;
			}
		}
	}

	private String[] carrierValues(int minCarrierFrequency, int maxCarrierFrequency) {
		int count = (maxCarrierFrequency - minCarrierFrequency) / 50 + 1;
		String[] values = new String[count];
		for (int i = 0; i < count; ++i)
			values[i] = String.format(Locale.US, "%d", i * 50 + minCarrierFrequency);
		return values;
	}

	private void setCarrierFrequency() {
		View view = getLayoutInflater().inflate(R.layout.carrier_frequency, null);
		NumberPicker picker = view.findViewById(R.id.carrier);
		int bandWidth = 1600;
		int maxCarrierFrequency = (outputRate - bandWidth) / 2;
		int minCarrierFrequency = outputChannel == 4 ? -maxCarrierFrequency : bandWidth / 2;
		if (carrierFrequency < minCarrierFrequency)
			carrierFrequency = minCarrierFrequency;
		if (carrierFrequency > maxCarrierFrequency)
			carrierFrequency = maxCarrierFrequency;
		picker.setMinValue(0);
		picker.setDisplayedValues(null);
		picker.setMaxValue((maxCarrierFrequency - minCarrierFrequency) / 50);
		picker.setValue((carrierFrequency - minCarrierFrequency) / 50);
		picker.setDisplayedValues(carrierValues(minCarrierFrequency, maxCarrierFrequency));
		setInputType(picker, InputType.TYPE_CLASS_NUMBER);
		AlertDialog.Builder builder = new AlertDialog.Builder(this, R.style.Theme_AlertDialog);
		builder.setTitle(R.string.carrier_frequency);
		builder.setView(view);
		builder.setNegativeButton(R.string.cancel, null);
		builder.setPositiveButton(R.string.okay, (dialogInterface, i) -> carrierFrequency = picker.getValue() * 50 + minCarrierFrequency);
		builder.show();
	}

	private void editCallSign() {
		View view = getLayoutInflater().inflate(R.layout.call_sign, null);
		EditText edit = view.findViewById(R.id.call);
		edit.setText(callSign);
		AlertDialog.Builder builder = new AlertDialog.Builder(this, R.style.Theme_AlertDialog);
		builder.setTitle(R.string.call_sign);
		builder.setView(view);
		builder.setNegativeButton(R.string.cancel, null);
		builder.setPositiveButton(R.string.okay, (dialogInterface, i) -> callSign = edit.getText().toString());
		builder.show();
	}

	private void showTextPage(String title, String message) {
		AlertDialog.Builder builder = new AlertDialog.Builder(this, R.style.Theme_AlertDialog);
		builder.setNeutralButton(R.string.close, null);
		builder.setTitle(title);
		builder.setMessage(message);
		builder.show();
	}

	@Override
	protected void onResume() {
		startListening();
		super.onResume();
	}

	@Override
	protected void onPause() {
		stopListening();
		storeSettings();
		super.onPause();
	}

	@Override
	protected void onDestroy() {
		audioTrack.stop();
		destroyEncoder();
		destroyDecoder();
		super.onDestroy();
	}
}