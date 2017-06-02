

import tensorflow as tf
import numpy as np


import cifar_input as cifar_inputs

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# GAN for cifar



image_size = 32
img_px = image_size*image_size
num_channels = 1

batch_size = 128

def plot(samples):
    fig = plt.figure(figsize=(4, 4))
    gs = gridspec.GridSpec(4, 4)
    gs.update(wspace=0.05, hspace=0.05)

    for i, sample in enumerate(samples):
        ax = plt.subplot(gs[i])
        plt.axis('off')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
        plt.imshow(sample.reshape(32, 32), cmap='Greys_r')

    return fig


def xavier_init(size):

    in_dim = size[0]
    xavier_stddev = 1./ tf.sqrt( in_dim/2.0 )
    return tf.random_normal( shape = size , stddev = xavier_stddev)


def build_model():
    # model requires discriminator and generator

    # discriminator parameters
    X = tf.placeholder( tf.float32 , shape=[None , img_px*num_channels ]   )
    # not random init of weights, but normal with some std
    D_w1 = tf.Variable( xavier_init( [ img_px*num_channels , batch_size ] )  )
    D_b1 = tf.Variable( tf.zeros( shape=[ batch_size ]) )

    D_w2 = tf.Variable( xavier_init( [ batch_size , 1 ] )  )
    D_b2 = tf.Variable( tf.zeros( shape=[ batch_size]) )

    theta_D = [ D_w1 , D_w2 , D_b1 , D_b2 ]


    # generator parameters

    Z = tf.placeholder( tf.float32 , shape =[ None , 100])
    G_w1 = tf.Variable( xavier_init( [ 100 , batch_size ] )  )
    G_b1 = tf.Variable( tf.zeros( shape=[ batch_size ]) )

    G_w2 = tf.Variable( xavier_init( [ batch_size , img_px*num_channels  ] )  )
    G_b2 = tf.Variable( tf.zeros( shape=[img_px*num_channels ]) )

    theta_G = [ G_w1 , G_w2 , G_b1 , G_b2 ]


    def sample_Z( m , n ):

        return np.random.uniform( -1.0 , 1 , size=[m ,n ])

    def generator(z):
        # generador 

        g_h1 = tf.nn.relu( tf.matmul( z , G_w1 ) +G_b1  )
        g_log_prob = tf.matmul( g_h1 , G_w2) + G_b2
        g_prob = tf.nn.sigmoid( g_log_prob )

        return g_prob

    def discriminator(x):

        d_h1 = tf.nn.relu( tf.matmul( x, D_w1) +  D_b1 )
        d_logit = tf.matmul( d_h1 , D_w2) + D_b2
        d_prob = tf.nn.sigmoid( d_logit )

        return d_prob , d_logit


    # build model 
    G_sample = generator(Z)

    D_real , D_logit_real = discriminator(X )
    D_fake , D_logit_fake = discriminator( G_sample )

    D_loss_real = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=D_logit_real, labels=tf.ones_like(D_logit_real)))
    D_loss_fake = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=D_logit_fake, labels=tf.zeros_like(D_logit_fake)))

    D_loss = D_loss_real + D_loss_fake
    G_loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=D_logit_fake, labels=tf.ones_like(D_logit_fake)))

    D_solver = tf.train.AdamOptimizer().minimize(D_loss , var_list=theta_D )
    G_solver = tf.train.AdamOptimizer().minimize( G_loss , var_list=theta_G )
    
    
    # train parameters 
    Z_dim = 100
    n_samples = 16 
    data_dir = "../data/cifar/cifar-10-py"

    images , cls , _ = cifar_inputs.load_training_data()

    images_1 = cifar_inputs.get_single_class( 5 , images , cls )

    batch = cifar_inputs.get_random_batch( images_1 ,  batch_size )
    
    #print( len(batch))
    
 
    #input handler 
    sess = tf.Session()
    sess.run(tf.global_variables_initializer() )
    steps = 150000
    for it in range( steps ):

        # generate n samples 
        #samples_g = sess.run( G_sample , feed_dict={ Z:sample_Z( n_samples , Z_dim) } )
        if it % 1000 == 0:
            print( "painting bitch" )
            samples_g = sess.run( G_sample , feed_dict={ Z:sample_Z( n_samples , Z_dim) } )
            fig = plot(samples_g)
            plt.savefig('./out/perro_{}.png'.format(it) , bbox_inches= 'tight' )
            plt.close(fig)
            
        print("#Training... ")
        #X_cifar , _ = cifar_inputs.inputs(False, data_dir, batch_size )
        X_np = cifar_inputs.get_random_batch(images_1 , batch_size ) 
        
        _ , D_loss_curr = sess.run( [ D_solver , D_loss ] , feed_dict={ X: X_np , Z : sample_Z( batch_size, Z_dim  ) }  ) 
        _ , G_loss_curr = sess.run( [G_solver, G_loss] , feed_dict={ Z: sample_Z( batch_size , Z_dim) } )
        
        
        if it % 1000 == 0:
            print('Iter: {}'.format(it))


            print('D loss: {:.4}'. format(D_loss_curr))
            print('G_loss: {:.4}'.format(G_loss_curr))
            print()

            

def main():

    build_model()

    print("aspero zocio")
    return

    
main()
