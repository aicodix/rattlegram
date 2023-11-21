package com.aicodix.rattlegram;

import android.Manifest;
import android.app.NotificationChannel;
import android.app.NotificationManager;
import android.app.PendingIntent;
import android.app.Service;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.os.IBinder;
import android.util.Log;

import androidx.annotation.Nullable;
import androidx.core.app.NotificationCompat;
import androidx.core.content.ContextCompat;

public class PassiveListenService extends Service {
    final String CHANNEL_ID = "Rattlegram Passive Listen";

    @Override
    public int onStartCommand(Intent intent, int flags, int startId) {
        int micPermission = ContextCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO);

        if (micPermission == PackageManager.PERMISSION_DENIED) {
            Log.e("PermissionError", "No permission to use mic");
            stopSelf();
        } else {
            NotificationChannel channel = new NotificationChannel(CHANNEL_ID, "Rattlegram", NotificationManager.IMPORTANCE_DEFAULT);
            channel.setDescription("Rattlegram channel for foreground service notifications");
            NotificationManager notificationManager = getSystemService(NotificationManager.class);
            notificationManager.createNotificationChannel(channel);

            Intent mainActIntent = new Intent(this, MainActivity.class);
            PendingIntent mainActPendingIntent = PendingIntent.getActivity(this, 0, mainActIntent, PendingIntent.FLAG_IMMUTABLE);
            NotificationCompat.Builder notification = new NotificationCompat.Builder(this, CHANNEL_ID)
                    .setContentText("Foreground service running")
                    .setContentTitle("Rattlegram Passive Listen")
                    .setSmallIcon(R.drawable.ic_launcher_foreground)
                    .setContentIntent(mainActPendingIntent);
            startForeground(1, notification.build());
        }
        return super.onStartCommand(intent, flags, startId);
    }

    @Nullable
    @Override
    public IBinder onBind(Intent intent) {
        return null;
    }
}
